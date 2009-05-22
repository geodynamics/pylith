// -*- C++ -*-
//
// ----------------------------------------------------------------------
//
//                           Brad T. Aagaard
//                        U.S. Geological Survey
//
// {LicenseText}
//
// ----------------------------------------------------------------------
//

#include <portinfo>

#include "DirichletBoundary.hh" // implementation of object methods

#include "pylith/topology/Mesh.hh" // USES Mesh
#include "pylith/topology/SubMesh.hh" // USES SubMesh
#include "pylith/topology/Field.hh" // USES Field
#include "pylith/topology/Fields.hh" // USES Fields
#include "spatialdata/geocoords/CoordSys.hh" // USES CoordSys
#include "spatialdata/units/Nondimensional.hh" // USES Nondimensional

#include <strings.h> // USES strcasecmp()
#include <cassert> // USES assert()
#include <stdexcept> // USES std::runtime_error
#include <sstream> // USES std::ostringstream

// ----------------------------------------------------------------------
// Default constructor.
pylith::bc::DirichletBoundary::DirichletBoundary(void) :
  _boundaryMesh(0),
  _fields(0)
{ // constructor
} // constructor

// ----------------------------------------------------------------------
// Destructor.
pylith::bc::DirichletBoundary::~DirichletBoundary(void)
{ // destructor
  delete _boundaryMesh; _boundaryMesh = 0;
  delete _fields; _fields = 0;
} // destructor

// ----------------------------------------------------------------------
// Initialize boundary condition.
void
pylith::bc::DirichletBoundary::initialize(const topology::Mesh& mesh,
					  const double upDir[3])
{ // initialize
  const int numFixedDOF = _fixedDOF.size();
  if (0 == numFixedDOF)
    return;

  _boundaryMesh = new topology::SubMesh(mesh, _label.c_str());
  assert(0 != _boundaryMesh);

  _getPoints(mesh);
  _setupQueryDatabases();
  _queryDatabases(mesh);
} // initialize

// ----------------------------------------------------------------------
// Get vertex field of BC initial or rate of change of values.
const pylith::topology::Field<pylith::topology::SubMesh>&
pylith::bc::DirichletBoundary::vertexField(const char* name,
					   const topology::SolutionFields& fields)
{ // getVertexField
  typedef topology::SubMesh::SieveMesh SieveMesh;
  typedef topology::SubMesh::RealSection RealSection;

  assert(0 != name);
  assert(0 != _boundaryMesh);
  assert(0 != _normalizer);

  const ALE::Obj<SieveMesh>& sieveMesh = _boundaryMesh->sieveMesh();
  assert(!sieveMesh.isNull());

  const ALE::Obj<SieveMesh::label_sequence>& vertices = 
    sieveMesh->depthStratum(0);
  assert(!vertices.isNull());

  const spatialdata::geocoords::CoordSys* cs = _boundaryMesh->coordsys();
  assert(0 != cs);
  const int fiberDim = cs->spaceDim();
  double_array values(fiberDim);

  const int numPoints = _points.size();
  const int numFixedDOF = _fixedDOF.size();

  if (0 == _fields) {
    _fields = 
      new topology::Fields<topology::Field<topology::SubMesh> >(*_boundaryMesh);
    assert(0 != _fields);
    _fields->add("output buffer", "temporary_buffer");
    topology::Field<topology::SubMesh>& buffer =
      _fields->get("output buffer");
    buffer.newSection(vertices, fiberDim);
    buffer.allocate();
    buffer.addDimensionOkay(true);
    buffer.vectorFieldType(topology::Field<topology::SubMesh>::VECTOR);
  } // if

  topology::Field<topology::SubMesh>& buffer =
    _fields->get("output buffer");
  if (0 == strcasecmp(name, "initial")) {
    buffer.label("bc_displacement");
    buffer.vectorFieldType(topology::Field<topology::SubMesh>::VECTOR);
    buffer.scale(_normalizer->lengthScale());
    buffer.zero();
    const ALE::Obj<RealSection>& section = buffer.section();

    for (int iPoint=0; iPoint < numPoints; ++iPoint) {
      const SieveMesh::point_type point = _points[iPoint];
      assert(fiberDim == section->getFiberDimension(point));
      for (int iDOF=0; iDOF < numFixedDOF; ++iDOF)
	values[_fixedDOF[iDOF]] = _valuesInitial[iPoint*numFixedDOF+iDOF];
      section->updatePointAll(_points[iPoint], &values[0]);
    } // for
  } else if (0 == strcasecmp(name, "rate-of-change")) {
    buffer.label("bc_velocity");
    buffer.vectorFieldType(topology::Field<topology::SubMesh>::VECTOR);
    assert(0 != _normalizer->timeScale());
    const double velocityScale =
      _normalizer->lengthScale() / _normalizer->timeScale();
    buffer.scale(velocityScale);
    buffer.zero();
    const ALE::Obj<RealSection>& section = buffer.section();

    for (int iPoint=0; iPoint < numPoints; ++iPoint) {
      const SieveMesh::point_type point = _points[iPoint];
      assert(fiberDim == section->getFiberDimension(point));
      for (int iDOF=0; iDOF < numFixedDOF; ++iDOF)
	values[_fixedDOF[iDOF]] = _valuesRate[iPoint*numFixedDOF+iDOF];
      section->updatePointAll(_points[iPoint], &values[0]);
    } // for
  } else {
    std::ostringstream msg;
    msg
      << "Unknown field '" << name << "' requested for Dirichlet boundary BC '" 
      << _label << "'.";
    throw std::runtime_error(msg.str());
  } // else

  return buffer;
} // getVertexField


// End of file 
