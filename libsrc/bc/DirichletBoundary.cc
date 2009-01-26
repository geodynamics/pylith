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

#include "pylith/topology/FieldUniform.hh" // USES FieldUniform
#include "pylith/topology/Mesh.hh" // USES Mesh
#include "spatialdata/geocoords/CoordSys.hh" // USES CoordSys
#include "spatialdata/units/Nondimensional.hh" // USES Nondimensional

#include <Selection.hh> // USES submesh algorithms

#include <strings.h> // USES strcasecmp()
#include <cassert> // USES assert()
#include <stdexcept> // USES std::runtime_error
#include <sstream> // USES std::ostringstream

// ----------------------------------------------------------------------
// Default constructor.
pylith::bc::DirichletBoundary::DirichletBoundary(void) :
  _tmpField(0)
{ // constructor
} // constructor

// ----------------------------------------------------------------------
// Destructor.
pylith::bc::DirichletBoundary::~DirichletBoundary(void)
{ // destructor
  delete _tmpField; _tmpField = 0;
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

  _createBoundaryMesh(mesh);
  _getPoints(mesh);
  _setupQueryDatabases();
  _queryDatabases(mesh);
} // initialize

// ----------------------------------------------------------------------
// Get vertex field of BC initial or rate of change of values.
const pylith::topology::Field&
pylith::bc::DirichletBoundary::vertexField(const char* name,
					   const topology::Mesh& mesh,
					   const topology::SolutionFields& fields)
{ // getVertexField
  assert(0 != name);
  assert(!_boundaryMesh.isNull());
  assert(0 != _normalizer);

  const ALE::Obj<SieveMesh::label_sequence>& vertices = 
    _boundaryMesh->depthStratum(0);
  assert(!vertices.isNull());
  const SieveMesh::label_sequence::iterator verticesEnd = vertices->end();

  const spatialdata::geocoords::CoordSys* cs = mesh.coordsys();
  assert(0 != cs);
  const int spaceDim = cs->spaceDim();
  double_array values(spaceDim);

  const int numPoints = _points.size();
  const int numFixedDOF = _fixedDOF.size();

  if (0 == _tmpField) {
    _tmpField = new topology::FieldUniform(_boundaryMesh, spaceDim);
    assert(0 != _tmpField);
    _tmpField->createSection(vertices);
  } // if

  if (0 == strcasecmp(name, "initial")) {
    _tmpField->name("displacement");
    _tmpField->vectorFieldType(topology::Field::VECTOR);
    _tmpField->scale(_normalizer->lengthScale());
    _tmpField->addDimensionOkay(true);
    _tmpField->zero();
    const ALE::Obj<SieveRealSection>& section = _tmpField->section();

    for (int iPoint=0; iPoint < numPoints; ++iPoint) {
      const SieveMesh::point_type point = _points[iPoint];
      assert(spaceDim == section->getFiberDimension(point));
      for (int iDOF=0; iDOF < numFixedDOF; ++iDOF)
	values[_fixedDOF[iDOF]] = _valuesInitial[iPoint*numFixedDOF+iDOF];
      section->updatePointAll(_points[iPoint], &values[0]);
    } // for
  } else if (0 == strcasecmp(name, "rate-of-change")) {
    _tmpField->name("velocity");
    _tmpField->vectorFieldType(topology::Field::VECTOR);
    _tmpField->scale(_normalizer->lengthScale());
    _tmpField->addDimensionOkay(true);
    _tmpField->zero();
    const ALE::Obj<SieveRealSection>& section = _tmpField->section();

    for (int iPoint=0; iPoint < numPoints; ++iPoint) {
      const SieveMesh::point_type point = _points[iPoint];
      assert(spaceDim == section->getFiberDimension(point));
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

  return *_tmpField;
} // getVertexField

// ----------------------------------------------------------------------
// Extract submesh associated with boundary.
void
pylith::bc::DirichletBoundary::_createBoundaryMesh(const topology::Mesh& mesh)
{ // _createBoundaryMesh
  const ALE::Obj<SieveMesh>& sieveMesh = mesh.sieveMesh();
  assert(!sieveMesh.isNull());

  const ALE::Obj<SieveMesh::int_section_type>& groupField = 
    sieveMesh->getIntSection(_label);
  if (groupField.isNull()) {
    std::ostringstream msg;
    msg << "Could not find group of points '" << _label << "' in mesh.";
    throw std::runtime_error(msg.str());
  } // if
  _boundaryMesh = 
    ALE::Selection<SieveMesh>::submeshV<SieveSubMesh>(sieveMesh, groupField);
  if (_boundaryMesh.isNull()) {
    std::ostringstream msg;
    msg << "Could not construct boundary mesh for Dirichlet boundary "
	<< "condition '" << _label << "'.";
    throw std::runtime_error(msg.str());
  } // if
  _boundaryMesh->setRealSection("coordinates", 
				sieveMesh->getRealSection("coordinates"));
  // Create the parallel overlap
  ALE::Obj<SieveSubMesh::send_overlap_type> sendParallelMeshOverlap =
    _boundaryMesh->getSendOverlap();
  ALE::Obj<SieveSubMesh::recv_overlap_type> recvParallelMeshOverlap =
    _boundaryMesh->getRecvOverlap();
  SieveMesh::renumbering_type& renumbering = sieveMesh->getRenumbering();
  //   Can I figure this out in a nicer way?
  ALE::SetFromMap<std::map<SieveMesh::point_type,SieveMesh::point_type> > globalPoints(renumbering);

  ALE::OverlapBuilder<>::constructOverlap(globalPoints, renumbering,
					  sendParallelMeshOverlap,
					  recvParallelMeshOverlap);
  _boundaryMesh->setCalculatedOverlap(true);
} // _createBoundaryMesh


// End of file 
