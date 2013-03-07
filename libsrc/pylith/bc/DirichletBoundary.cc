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
// Copyright (c) 2010-2012 University of California, Davis
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
					  const PylithScalar upDir[3])
{ // initialize
  DirichletBC::initialize(mesh, upDir);

  _boundaryMesh = new topology::SubMesh(mesh, _label.c_str());
  assert(_boundaryMesh);
} // initialize

// ----------------------------------------------------------------------
// Get vertex field of BC initial or rate of change of values.
const pylith::topology::Field<pylith::topology::SubMesh>&
pylith::bc::DirichletBoundary::vertexField(const char* name,
					   const topology::SolutionFields& fields)
{ // getVertexField
  assert(_normalizer);
  const PylithScalar lengthScale = _normalizer->lengthScale();
  const PylithScalar timeScale = _normalizer->timeScale();
  const PylithScalar rateScale = lengthScale / timeScale;

  assert(_boundaryMesh);
  const spatialdata::geocoords::CoordSys* cs = _boundaryMesh->coordsys();
  assert(cs);
  const int spaceDim = cs->spaceDim();

  ALE::MemoryLogger& logger = ALE::MemoryLogger::singleton();
  logger.stagePush("BoundaryConditions");

  if (0 == _outputFields)
    _outputFields = new topology::Fields<topology::Field<topology::SubMesh> >(*_boundaryMesh);
  assert(_outputFields);
  _outputFields->add("buffer (vector)", "buffer_vector", topology::FieldBase::CELLS_FIELD, spaceDim);
  _outputFields->get("buffer (vector)").vectorFieldType(topology::FieldBase::VECTOR);
  _outputFields->get("buffer (vector)").scale(lengthScale);
  _outputFields->get("buffer (vector)").allocate();
  _outputFields->add("buffer (scalar)", "buffer_scalar", topology::FieldBase::CELLS_FIELD, 1);
  _outputFields->get("buffer (scalar)").vectorFieldType(topology::FieldBase::SCALAR);
  _outputFields->get("buffer (scalar)").scale(timeScale);
  _outputFields->get("buffer (scalar)").allocate();

  logger.stagePop();

  if (0 == strcasecmp(name, "initial_value"))
    return _bufferVector("initial", "initial_displacement", lengthScale);
  else if (0 == strcasecmp(name, "rate_of_change"))
    return _bufferVector("rate", "velocity", rateScale);
  else if (0 == strcasecmp(name, "change_in_value"))
    return _bufferVector("change", "displacement_change", lengthScale);
  else if (0 == strcasecmp(name, "rate_start_time"))
    return _bufferScalar("rate time", "velocity_start_time", timeScale);
  else if (0 == strcasecmp(name, "change_start_time"))
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
					     const PylithScalar scale)
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
  
  const spatialdata::geocoords::CoordSys* cs = _boundaryMesh->coordsys();
  assert(cs);
  const int spaceDim = cs->spaceDim();

  const int numPoints = _points.size();
  const int numFixedDOF = _bcDOF.size();
  PetscErrorCode err = 0;

  assert(_outputFields->hasField("buffer (vector)"));
  PetscSection outputSection = _outputFields->get("buffer (vector)").petscSection();assert(outputSection);
  PetscVec outputVec = _outputFields->get("buffer (vector)").localVector();assert(outputVec);
  PetscScalar *outputArray = NULL;
  err = VecGetArray(outputVec, &outputArray);CHECK_PETSC_ERROR(err);

  PetscSection fieldSection = _parameters->get(name).petscSection();assert(fieldSection);
  PetscVec fieldVec = _parameters->get(name).localVector();assert(fieldVec);
  PetscScalar *fieldArray = NULL;
  err = VecGetArray(fieldVec, &fieldArray);CHECK_PETSC_ERROR(err);

  for (int iPoint=0; iPoint < numPoints; ++iPoint) {
    const SieveMesh::point_type point = _points[iPoint];

    PetscInt odof, ooff, fdof, foff;
    err = PetscSectionGetDof(outputSection, point, &odof);CHECK_PETSC_ERROR(err);
    err = PetscSectionGetOffset(outputSection, point, &ooff);CHECK_PETSC_ERROR(err);
    err = PetscSectionGetDof(fieldSection, point, &fdof);CHECK_PETSC_ERROR(err);
    err = PetscSectionGetOffset(fieldSection, point, &foff);CHECK_PETSC_ERROR(err);
    assert(spaceDim == odof);
    assert(fdof == numFixedDOF);

    for(PetscInt iDOF=0; iDOF < numFixedDOF; ++iDOF)
      outputArray[ooff+_bcDOF[iDOF]] = fieldArray[foff+iDOF];
  } // for
  err = VecRestoreArray(outputVec, &outputArray);CHECK_PETSC_ERROR(err);
  err = VecRestoreArray(fieldVec, &fieldArray);CHECK_PETSC_ERROR(err);

  topology::Field<topology::SubMesh>& buffer = _outputFields->get("buffer (vector)");  
  buffer.label(label);
  buffer.scale(scale);

  return buffer;
} // _bufferVector

// ----------------------------------------------------------------------
// Get vertex scalar field with BC information.
const pylith::topology::Field<pylith::topology::SubMesh>&
pylith::bc::DirichletBoundary::_bufferScalar(const char* name,
					     const char* label,
					     const PylithScalar scale)
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
  
  const int numPoints = _points.size();
  PetscErrorCode err = 0;

  assert(_outputFields->hasField("buffer (scalar)"));
  PetscSection outputSection = _outputFields->get("buffer (scalar)").petscSection();assert(outputSection);
  PetscVec outputVec = _outputFields->get("buffer (scalar)").localVector();assert(outputVec);
  PetscScalar *outputArray = NULL;
  err = VecGetArray(outputVec, &outputArray);CHECK_PETSC_ERROR(err);
  
  PetscSection fieldSection = _parameters->get(name).petscSection();assert(fieldSection);
  PetscVec fieldVec = _parameters->get(name).localVector();assert(fieldVec);
  PetscScalar *fieldArray = NULL;
  err = VecGetArray(fieldVec, &fieldArray);CHECK_PETSC_ERROR(err);

  for (int iPoint=0; iPoint < numPoints; ++iPoint) {
    const SieveMesh::point_type point = _points[iPoint];

    PetscInt odof, ooff, fdof, foff;
    err = PetscSectionGetDof(outputSection, point, &odof);CHECK_PETSC_ERROR(err);
    err = PetscSectionGetOffset(outputSection, point, &ooff);CHECK_PETSC_ERROR(err);
    err = PetscSectionGetDof(fieldSection, point, &fdof);CHECK_PETSC_ERROR(err);
    err = PetscSectionGetOffset(fieldSection, point, &foff);CHECK_PETSC_ERROR(err);
    assert(1 == odof);
    assert(fdof == 1);

    outputArray[ooff] = fieldArray[foff];
  } // for
  
  topology::Field<topology::SubMesh>& buffer = _outputFields->get("buffer (scalar)");  
  buffer.label(label);
  buffer.scale(scale);

  return buffer;
} // _bufferScalar


// End of file 
