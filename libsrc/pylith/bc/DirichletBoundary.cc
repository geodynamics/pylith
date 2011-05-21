// -*- C++ -*-
//
// ----------------------------------------------------------------------
//
// Brad T. Aagaard, U.S. Geological Survey
// Charles A. Williams, GNS Science
// Matthew G. Knepley, University of Chicago
//
// This code was developed as part of the Computational Infrastructure
// for Geodynamics (http://geodynamics.org).
//
// Copyright (c) 2010 University of California, Davis
//
// See COPYING for license information.
//
// ----------------------------------------------------------------------
//

#include <portinfo>

#include "DirichletBoundary.hh" // implementation of object methods

#include "pylith/topology/Mesh.hh" // USES Mesh
#include "pylith/topology/SubMesh.hh" // USES SubMesh
#include "pylith/topology/Field.hh" // USES Field
#include "pylith/topology/FieldsNew.hh" // USES FieldsNew
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
  deallocate();
} // destructor

// ----------------------------------------------------------------------
// Deallocate PETSc and local data structures.
void
pylith::bc::DirichletBoundary::deallocate(void)
{ // deallocate
  DirichletBC::deallocate();

  delete _boundaryMesh; _boundaryMesh = 0;
  delete _outputFields; _outputFields = 0;
} // deallocate
  
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
  assert(_normalizer);
  const double lengthScale = _normalizer->lengthScale();
  const double timeScale = _normalizer->timeScale();
  const double rateScale = lengthScale / timeScale;

  assert(_boundaryMesh);
  const spatialdata::geocoords::CoordSys* cs = _boundaryMesh->coordsys();
  assert(0 != cs);
  const int spaceDim = cs->spaceDim();

  ALE::MemoryLogger& logger = ALE::MemoryLogger::singleton();
  logger.stagePush("BoundaryConditions");

  if (0 == _outputFields)
    _outputFields = new topology::FieldsNew<topology::SubMesh>(*_boundaryMesh);
  assert(0 != _outputFields);
  _outputFields->add("buffer (vector)", "buffer_vector", 
		     spaceDim, 
		     topology::FieldBase::VECTOR,
		     lengthScale);
  _outputFields->add("buffer (scalar)", "buffer_scalar", 
		     1, 
		     topology::FieldBase::SCALAR,
		     timeScale);
  _outputFields->allocate(topology::FieldBase::CELLS_FIELD, 1);

  logger.stagePop();

  if (0 == strcasecmp(name, "initial-value"))
    return _bufferVector("initial", "initial_displacement", lengthScale);
  else if (0 == strcasecmp(name, "rate-of-change"))
    return _bufferVector("rate", "velocity", rateScale);
  else if (0 == strcasecmp(name, "change-in-value"))
    return _bufferVector("change", "displacement_change", lengthScale);
  else if (0 == strcasecmp(name, "rate-start-time"))
    return _bufferScalar("rate time", "velocity_start_time", timeScale);
  else if (0 == strcasecmp(name, "change-start-time"))
    return _bufferScalar("change time", "change_start_time", timeScale);
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
  typedef topology::Mesh::RealUniformSection RealUniformSection;
  typedef topology::SubMesh::RealSection SubRealSection;
  typedef topology::SubMesh::RealUniformSection SubRealUniformSection;

  assert(_boundaryMesh);
  assert(_parameters);
  assert(_outputFields);

  if (!_parameters->hasField(name)) {
    std::ostringstream msg;
    msg << "Parameters for field '" << label << " not provided in "
	<< "Dirichlet BC '" << _label << "'.";
    throw std::runtime_error(msg.str());
  } // if
  
  const ALE::Obj<SieveMesh>& sieveMesh = _boundaryMesh->sieveMesh();
  assert(!sieveMesh.isNull());

  const spatialdata::geocoords::CoordSys* cs = _boundaryMesh->coordsys();
  assert(0 != cs);
  const int spaceDim = cs->spaceDim();

  const int numPoints = _points.size();
  const int numFixedDOF = _bcDOF.size();

  assert(_outputFields->hasField("buffer (vector)"));
  const ALE::Obj<SubRealUniformSection>& outputSection = 
    _outputFields->section();
  assert(!outputSection.isNull());
  const int outputFiberDim = _outputFields->fiberDim();  
  double_array outputVertex(outputFiberDim);
  const int bufferIndex = _outputFields->sectionIndex("buffer (vector)");
  const int bufferFiberDim = _outputFields->sectionFiberDim("buffer (vector)");
  assert(bufferIndex + bufferFiberDim <= outputFiberDim);
  assert(spaceDim == bufferFiberDim);
  
  const ALE::Obj<RealUniformSection>& parametersSection = 
    _parameters->section();
  assert(!parametersSection.isNull());
  const int parametersFiberDim = _parameters->fiberDim();
  const int fieldIndex = _parameters->sectionIndex(name);
  const int fieldFiberDim = _parameters->sectionFiberDim(name);
  assert(fieldIndex + fieldFiberDim <= parametersFiberDim);
  assert(fieldFiberDim == numFixedDOF);
  
  for (int iPoint=0; iPoint < numPoints; ++iPoint) {
    const SieveMesh::point_type point = _points[iPoint];
    outputVertex = 0.0;

    assert(parametersFiberDim == parametersSection->getFiberDimension(point));
    const double* parametersVertex = parametersSection->restrictPoint(point);
    assert(parametersVertex);
    
    for (int iDOF=0; iDOF < numFixedDOF; ++iDOF)
      outputVertex[bufferIndex+_bcDOF[iDOF]] = 
	parametersVertex[fieldIndex+iDOF];
    assert(outputFiberDim == outputSection->getFiberDimension(point));
    outputSection->updatePointAll(point, &outputVertex[0]);
  } // for

  topology::Field<topology::SubMesh>& buffer =
    _outputFields->get("buffer (vector)");  
  buffer.label(label);
  buffer.scale(scale);

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
  typedef topology::Mesh::RealUniformSection RealUniformSection;
  typedef topology::SubMesh::RealSection SubRealSection;
  typedef topology::SubMesh::RealUniformSection SubRealUniformSection;

  assert(_boundaryMesh);
  assert(_parameters);
  assert(_outputFields);

  if (!_parameters->hasField(name)) {
    std::ostringstream msg;
    msg << "Parameters for field '" << label << " not provided in "
	<< "Dirichlet BC '" << _label << "'.";
    throw std::runtime_error(msg.str());
  } // if
  
  const ALE::Obj<SieveMesh>& sieveMesh = _boundaryMesh->sieveMesh();
  assert(!sieveMesh.isNull());

  const int numPoints = _points.size();

  assert(_outputFields->hasField("buffer (scalar)"));
  const ALE::Obj<SubRealUniformSection>& outputSection = 
    _outputFields->section();
  assert(!outputSection.isNull());
  const int outputFiberDim = _outputFields->fiberDim();  
  double_array outputVertex(outputFiberDim);
  const int bufferIndex = _outputFields->sectionIndex("buffer (vector)");
  const int bufferFiberDim = _outputFields->sectionFiberDim("buffer (vector)");
  assert(bufferIndex + bufferFiberDim <= outputFiberDim);
  assert(1 == bufferFiberDim);
  
  const ALE::Obj<RealUniformSection>& parametersSection = 
    _parameters->section();
  assert(!parametersSection.isNull());
  const int parametersFiberDim = _parameters->fiberDim();
  const int fieldIndex = _parameters->sectionIndex(name);
  const int fieldFiberDim = _parameters->sectionFiberDim(name);
  assert(fieldIndex + fieldFiberDim <= parametersFiberDim);
  assert(1 == fieldFiberDim);
  
  for (int iPoint=0; iPoint < numPoints; ++iPoint) {
    const SieveMesh::point_type point = _points[iPoint];
    outputVertex = 0.0;

    assert(parametersFiberDim == parametersSection->getFiberDimension(point));
    const double* parametersVertex = parametersSection->restrictPoint(point);
    assert(parametersVertex);
    
    outputVertex[bufferIndex] = parametersVertex[fieldIndex];
    assert(outputFiberDim == outputSection->getFiberDimension(point));
    outputSection->updatePointAll(point, &outputVertex[0]);
  } // for
  
  topology::Field<topology::SubMesh>& buffer =
    _outputFields->get("buffer (scalar)");  
  buffer.label(label);
  buffer.scale(scale);

  return buffer;
} // _bufferScalar


// End of file 
