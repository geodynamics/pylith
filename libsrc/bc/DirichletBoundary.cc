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
  _outputFields(0)
{ // constructor
} // constructor

// ----------------------------------------------------------------------
// Destructor.
pylith::bc::DirichletBoundary::~DirichletBoundary(void)
{ // destructor
  delete _boundaryMesh; _boundaryMesh = 0;
  delete _outputFields; _outputFields = 0;
} // destructor

// ----------------------------------------------------------------------
// Initialize boundary condition.
void
pylith::bc::DirichletBoundary::initialize(const topology::Mesh& mesh,
					  const double upDir[3])
{ // initialize
  DirichletBC::initialize(mesh, upDir);

  _boundaryMesh = new topology::SubMesh(mesh, _label.c_str());
  assert(0 != _boundaryMesh);
} // initialize

// ----------------------------------------------------------------------
// Get vertex field of BC initial or rate of change of values.
const pylith::topology::Field<pylith::topology::SubMesh>&
pylith::bc::DirichletBoundary::vertexField(const char* name,
					   const topology::SolutionFields& fields)
{ // getVertexField
  assert(0 != _normalizer);
  const double lengthScale = _normalizer->lengthScale();
  const double timeScale = _normalizer->timeScale();
  const double rateScale = lengthScale / timeScale;

  if (0 == strcasecmp(name, "initial"))
    return _bufferVector("initial", "initial_displacement", lengthScale);
  else if (0 == strcasecmp(name, "rate-of-change"))
    return _bufferVector("rate", "initial_velocity", rateScale);
  else if (0 == strcasecmp(name, "change-in-value"))
    return _bufferVector("change", "displacement_change", lengthScale);
  else if (0 == strcasecmp(name, "rate-start-time"))
    return _bufferScalar("rate-start-time", "velocity_start_time", timeScale);
  else if (0 == strcasecmp(name, "change-start-time"))
    return _bufferScalar("change-start-time", "change_start_time", timeScale);
  else {
    std::ostringstream msg;
    msg
      << "Unknown field '" << name << "' requested for Dirichlet boundary BC '" 
      << _label << "'.";
    throw std::runtime_error(msg.str());
  } // else

  // Satisfy return value (should never reach here)
  return _outputFields->get("null");
} // getVertexField

// ----------------------------------------------------------------------
// Get vertex vector field with BC information.
const pylith::topology::Field<pylith::topology::SubMesh>&
pylith::bc::DirichletBoundary::_bufferVector(const char* name,
					     const char* label,
					     const double scale)
{ // _bufferVector
  typedef topology::SubMesh::SieveMesh SieveMesh;
  typedef topology::SubMesh::RealSection RealSection;

  assert(0 != _boundaryMesh);
  assert(0 != _parameters);

  if (_parameters->hasField(name)) {
    std::ostringstream msg;
    msg << "Parameters for field '" << label << " not provided in "
	<< "Dirichlet BC '" << _label << "'.";
    throw std::runtime_error(msg.str());
  } // if
  
  if (0 == _outputFields)
    _outputFields = 
      new topology::Fields<topology::Field<topology::SubMesh> >(*_boundaryMesh);
  assert(0 != _outputFields);
  
  const ALE::Obj<SieveMesh>& sieveMesh = _boundaryMesh->sieveMesh();
  assert(!sieveMesh.isNull());

  const spatialdata::geocoords::CoordSys* cs = _boundaryMesh->coordsys();
  assert(0 != cs);
  const int fiberDim = cs->spaceDim();

  const int numPoints = _points.size();
  const int numFixedDOF = _bcDOF.size();
  ALE::MemoryLogger& logger = ALE::MemoryLogger::singleton();
  logger.stagePush("BoundaryConditions");

  double_array bufferVertex(fiberDim);
  if (!_outputFields->hasField("buffer (vector)")) {
    _outputFields->add("buffer (vector)", "buffer");
    topology::Field<topology::SubMesh>& buffer =
      _outputFields->get("buffer (vector)");  
    buffer.newSection(topology::FieldBase::VERTICES_FIELD, fiberDim);
    buffer.allocate();
  } // if
  topology::Field<topology::SubMesh>& buffer =
    _outputFields->get("buffer (vector)");  
  buffer.label(label);
  buffer.scale(scale);
  buffer.vectorFieldType(topology::FieldBase::VECTOR);
  buffer.zero();
  const ALE::Obj<RealSection>& bufferSection = buffer.section();
  assert(!bufferSection.isNull());

  double_array parameterVertex(numFixedDOF);
  const ALE::Obj<RealSection>& parameterSection = 
    _parameters->get(name).section();
  assert(!parameterSection.isNull());
  
  for (int iPoint=0; iPoint < numPoints; ++iPoint) {
    const SieveMesh::point_type point = _points[iPoint];
    bufferVertex = 0.0;

    parameterSection->restrictPoint(point, &parameterVertex[0], 
				    parameterVertex.size());
    
    for (int iDOF=0; iDOF < numFixedDOF; ++iDOF)
      bufferVertex[_bcDOF[iDOF]] = parameterVertex[iDOF];
    assert(fiberDim == bufferSection->getFiberDimension(point));
    bufferSection->updatePointAll(point, &bufferVertex[0]);
  } // for

  logger.stagePop();
  return buffer;
} // _bufferVector

// ----------------------------------------------------------------------
// Get vertex scalar field with BC information.
const pylith::topology::Field<pylith::topology::SubMesh>&
pylith::bc::DirichletBoundary::_bufferScalar(const char* name,
					     const char* label,
					     const double scale)
{ // _bufferScalar
  typedef topology::SubMesh::SieveMesh SieveMesh;
  typedef topology::SubMesh::RealSection RealSection;

  assert(0 != _boundaryMesh);
  assert(0 != _parameters);

  if (_parameters->hasField(name)) {
    std::ostringstream msg;
    msg << "Parameters for field '" << label << " not provided in "
	<< "Dirichlet BC '" << _label << "'.";
    throw std::runtime_error(msg.str());
  } // if
  
  if (0 == _outputFields)
    _outputFields = 
      new topology::Fields<topology::Field<topology::SubMesh> >(*_boundaryMesh);
  assert(0 != _outputFields);
  
  const ALE::Obj<SieveMesh>& sieveMesh = _boundaryMesh->sieveMesh();
  assert(!sieveMesh.isNull());

  const int numPoints = _points.size();
  const int fiberDim = 1;

  if (!_outputFields->hasField("buffer (scalar)")) {
    _outputFields->add("buffer (scalar)", "buffer");
    topology::Field<topology::SubMesh>& buffer =
      _outputFields->get("buffer (scalar)");  
    buffer.newSection(topology::FieldBase::VERTICES_FIELD, fiberDim);
    buffer.allocate();
  } // if
  topology::Field<topology::SubMesh>& buffer =
    _outputFields->get("buffer (scalar)");  
  buffer.label(label);
  buffer.scale(scale);
  buffer.vectorFieldType(topology::FieldBase::SCALAR);
  buffer.zero();
  const ALE::Obj<RealSection>& bufferSection = buffer.section();
  assert(!bufferSection.isNull());

  const ALE::Obj<RealSection>& parameterSection = 
    _parameters->get(name).section();
  assert(!parameterSection.isNull());
  
  for (int iPoint=0; iPoint < numPoints; ++iPoint) {
    const SieveMesh::point_type point = _points[iPoint];

    assert(1 == bufferSection->getFiberDimension(point));
    assert(1 == parameterSection->getFiberDimension(point));
    bufferSection->updatePointAll(point, 
				  parameterSection->restrictPoint(point));
  } // for
  
  return buffer;
} // _bufferScalar


// End of file 
