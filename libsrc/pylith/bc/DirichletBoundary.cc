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
// Copyright (c) 2010-2014 University of California, Davis
//
// See COPYING for license information.
//
// ----------------------------------------------------------------------
//

#include <portinfo>

#include "DirichletBoundary.hh" // implementation of object methods

#include "pylith/topology/Mesh.hh" // USES Mesh
#include "pylith/topology/Field.hh" // USES Field
#include "pylith/topology/Fields.hh" // USES Fields
#include "pylith/topology/VisitorMesh.hh" // USES VecVisitorMesh
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
  PYLITH_METHOD_BEGIN;

  DirichletBC::deallocate();

  delete _boundaryMesh; _boundaryMesh = 0;
  delete _outputFields; _outputFields = 0;

  PYLITH_METHOD_END;
} // deallocate
  
// ----------------------------------------------------------------------
// Initialize boundary condition.
void
pylith::bc::DirichletBoundary::initialize(const topology::Mesh& mesh,
					  const PylithScalar upDir[3])
{ // initialize
  PYLITH_METHOD_BEGIN;

  DirichletBC::initialize(mesh, upDir);

  _boundaryMesh = new topology::Mesh(mesh, _label.c_str());
  assert(_boundaryMesh);

  PYLITH_METHOD_END;
} // initialize

// ----------------------------------------------------------------------
// Get vertex field of BC initial or rate of change of values.
const pylith::topology::Field&
pylith::bc::DirichletBoundary::vertexField(const char* name,
					   const topology::SolutionFields& fields)
{ // getVertexField
  PYLITH_METHOD_BEGIN;

  assert(_normalizer);
  const PylithScalar lengthScale = _normalizer->lengthScale();
  const PylithScalar timeScale = _normalizer->timeScale();
  const PylithScalar rateScale = lengthScale / timeScale;

  assert(_boundaryMesh);
  const spatialdata::geocoords::CoordSys* cs = _boundaryMesh->coordsys();
  assert(cs);
  const int spaceDim = cs->spaceDim();

  if (!_outputFields) {
    _outputFields = new topology::Fields(*_boundaryMesh);
  } // if
  assert(_outputFields);
  _outputFields->add("buffer (vector)", "buffer_vector", topology::FieldBase::FACES_FIELD, spaceDim);
  topology::Field& bufferVector = _outputFields->get("buffer (vector)");
  bufferVector.vectorFieldType(topology::FieldBase::VECTOR);
  bufferVector.scale(lengthScale);
  bufferVector.allocate();

  _outputFields->add("buffer (scalar)", "buffer_scalar", topology::FieldBase::FACES_FIELD, 1);
  topology::Field& bufferScalar = _outputFields->get("buffer (scalar)");
  bufferScalar.vectorFieldType(topology::FieldBase::SCALAR);
  bufferScalar.scale(timeScale);
  bufferScalar.allocate();

  if (0 == strcasecmp(name, "initial_value"))
    PYLITH_METHOD_RETURN(_bufferVector("initial", "initial_displacement", lengthScale));
  else if (0 == strcasecmp(name, "rate_of_change"))
    PYLITH_METHOD_RETURN(_bufferVector("rate", "velocity", rateScale));
  else if (0 == strcasecmp(name, "change_in_value"))
    PYLITH_METHOD_RETURN(_bufferVector("change", "displacement_change", lengthScale));
  else if (0 == strcasecmp(name, "rate_start_time"))
    PYLITH_METHOD_RETURN(_bufferScalar("rate time", "velocity_start_time", timeScale));
  else if (0 == strcasecmp(name, "change_start_time"))
    PYLITH_METHOD_RETURN(_bufferScalar("change time", "change_start_time", timeScale));
  else {
    std::ostringstream msg;
    msg
      << "Unknown field '" << name << "' requested for Dirichlet boundary BC '" 
      << _label << "'.";
    throw std::runtime_error(msg.str());
  } // else

  // Satisfy return value (should never reach here)

  PYLITH_METHOD_RETURN(_outputFields->get("null"));
} // getVertexField

// ----------------------------------------------------------------------
// Get vertex vector field with BC information.
const pylith::topology::Field&
pylith::bc::DirichletBoundary::_bufferVector(const char* name,
					     const char* label,
					     const PylithScalar scale)
{ // _bufferVector
  PYLITH_METHOD_BEGIN;

  assert(_boundaryMesh);
  assert(_parameters);
  assert(_outputFields);

  if (!_parameters->hasField(name)) {
    std::ostringstream msg;
    msg << "Parameters for field '" << label << " not provided in "
	<< "Dirichlet BC '" << _label << "'.";
    throw std::runtime_error(msg.str());
  } // if
  
  assert(_outputFields->hasField("buffer (vector)"));
  topology::Field& buffer = _outputFields->get("buffer (vector)");
  topology::VecVisitorMesh bufferVisitor(buffer);
  PetscScalar* bufferArray = bufferVisitor.localArray();

  topology::Field& field = _parameters->get(name);
  topology::VecVisitorMesh fieldVisitor(field);
  PetscScalar* fieldArray = fieldVisitor.localArray();

  const spatialdata::geocoords::CoordSys* cs = _boundaryMesh->coordsys();
  assert(cs);
  const int spaceDim = cs->spaceDim();

  const int numPoints = _points.size();
  const int numFixedDOF = _bcDOF.size();
  for (int iPoint=0; iPoint < numPoints; ++iPoint) {
    const PetscInt point = _points[iPoint];

    const PetscInt boff = bufferVisitor.sectionOffset(point);
    const PetscInt foff = fieldVisitor.sectionOffset(point);
    assert(spaceDim == bufferVisitor.sectionDof(point));
    assert(numFixedDOF == fieldVisitor.sectionDof(point));

    for(PetscInt iDOF=0; iDOF < numFixedDOF; ++iDOF)
      bufferArray[boff+_bcDOF[iDOF]] = fieldArray[foff+iDOF];
  } // for

  buffer.label(label);
  buffer.scale(scale);

  PYLITH_METHOD_RETURN(buffer);
} // _bufferVector

// ----------------------------------------------------------------------
// Get vertex scalar field with BC information.
const pylith::topology::Field&
pylith::bc::DirichletBoundary::_bufferScalar(const char* name,
					     const char* label,
					     const PylithScalar scale)
{ // _bufferScalar
  PYLITH_METHOD_BEGIN;

  assert(_boundaryMesh);
  assert(_parameters);
  assert(_outputFields);

  if (!_parameters->hasField(name)) {
    std::ostringstream msg;
    msg << "Parameters for field '" << label << "' not provided in "
	<< "Dirichlet BC '" << _label << "'.";
    throw std::runtime_error(msg.str());
  } // if
  
  assert(_outputFields->hasField("buffer (scalar)"));
  topology::Field& buffer = _outputFields->get("buffer (scalar)");
  topology::VecVisitorMesh bufferVisitor(buffer);
  PetscScalar* bufferArray = bufferVisitor.localArray();

  topology::Field& field = _parameters->get(name);
  topology::VecVisitorMesh fieldVisitor(field);
  PetscScalar* fieldArray = fieldVisitor.localArray();

  const int numPoints = _points.size();
  for (int iPoint=0; iPoint < numPoints; ++iPoint) {
    const PetscInt point = _points[iPoint];

    const PetscInt boff = bufferVisitor.sectionOffset(point);
    const PetscInt foff = fieldVisitor.sectionOffset(point);
    assert(1 == bufferVisitor.sectionDof(point));
    assert(1 == fieldVisitor.sectionDof(point));

    bufferArray[boff] = fieldArray[foff];
  } // for
  
  buffer.label(label);
  buffer.scale(scale);

  PYLITH_METHOD_RETURN(buffer);
} // _bufferScalar


// End of file 
