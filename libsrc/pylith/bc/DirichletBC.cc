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

#include "DirichletBC.hh" // implementation of object methods

#include "pylith/topology/Mesh.hh" // USES Mesh
#include "pylith/topology/Field.hh" // USES Field
#include "pylith/topology/Fields.hh" // USES Fields
#include "pylith/topology/VisitorMesh.hh" // USES VecVisitorMesh
#include "spatialdata/geocoords/CoordSys.hh" // USES CoordSys
#include "spatialdata/units/Nondimensional.hh" // USES Nondimensional

#include <cstring> // USES strcpy()
#include <cassert> // USES assert()
#include <stdexcept> // USES std::runtime_error
#include <sstream> // USES std::ostringstream

// ----------------------------------------------------------------------
// Default constructor.
pylith::bc::DirichletBC::DirichletBC(void)
{ // constructor
} // constructor

// ----------------------------------------------------------------------
// Destructor.
pylith::bc::DirichletBC::~DirichletBC(void)
{ // destructor
  deallocate();
} // destructor

// ----------------------------------------------------------------------
// Deallocate PETSc and local data structures.
void
pylith::bc::DirichletBC::deallocate(void)
{ // deallocate
  PYLITH_METHOD_BEGIN;

  TimeDependentPoints::deallocate();
  feassemble::Constraint::deallocate();

  PYLITH_METHOD_END;
} // deallocate
  
// ----------------------------------------------------------------------
// Initialize boundary condition.
void
pylith::bc::DirichletBC::initialize(const topology::Mesh& mesh,
				    const PylithScalar upDir[3])
{ // initialize
  PYLITH_METHOD_BEGIN;

  if (0 == _bcDOF.size())
    PYLITH_METHOD_END;

  _getPoints(mesh);

  assert(_normalizer);
  const PylithScalar lengthScale = _normalizer->lengthScale();
  _queryDatabases(mesh, lengthScale, "displacement");

  PYLITH_METHOD_END;
} // initialize

// ----------------------------------------------------------------------
// Set number of degrees of freedom that are constrained at points in field.
void
pylith::bc::DirichletBC::setConstraintSizes(const topology::Field& field)
{ // setConstraintSizes
  PYLITH_METHOD_BEGIN;

  const int numFixedDOF = _bcDOF.size();
  if (0 == numFixedDOF)
    PYLITH_METHOD_END;

  PetscSection section = field.localSection();assert(section);
  PetscInt numFields;
  PetscErrorCode err = 0;

  // Set constraints in field
  const int numPoints = _points.size();
  _offsetLocal.resize(numPoints);
  err = PetscSectionGetNumFields(section, &numFields);PYLITH_CHECK_ERROR(err);

  for (int iPoint=0; iPoint < numPoints; ++iPoint) {
    const PetscInt point = _points[iPoint];

    PetscInt dof, cdof;
    err = PetscSectionGetDof(section, point, &dof);PYLITH_CHECK_ERROR(err);
    err = PetscSectionGetConstraintDof(section, point, &cdof);PYLITH_CHECK_ERROR(err);
    if (cdof + numFixedDOF > dof) {
      std::ostringstream msg;
      msg
	<< "Found overly constrained point while setting up constraints for "
	<< "DirichletBC boundary condition '" << _label << "'. "
	<< "Number of DOF at point " << point << " is " << dof
	<< " and number of attempted constraints is " << cdof+numFixedDOF << ".";
      throw std::runtime_error(msg.str());
    } // if
    _offsetLocal[iPoint] = cdof;
    err = PetscSectionAddConstraintDof(section, point, numFixedDOF);PYLITH_CHECK_ERROR(err);
    // We should be specifying what field the BC is for
    if (numFields) {err = PetscSectionAddFieldConstraintDof(section, point, 0, numFixedDOF);PYLITH_CHECK_ERROR(err);}
  } // for

  PYLITH_METHOD_END;
} // setConstraintSizes

// ----------------------------------------------------------------------
// Set which degrees of freedom are constrained at points in field.
void
pylith::bc::DirichletBC::setConstraints(const topology::Field& field)
{ // setConstraints
  PYLITH_METHOD_BEGIN;

  const int numFixedDOF = _bcDOF.size();
  if (0 == numFixedDOF)
    PYLITH_METHOD_END;

  PetscSection section = field.localSection();assert(section);
  PetscInt numFields;
  PetscErrorCode err = 0;

  const int numPoints = _points.size();
  err = PetscSectionGetNumFields(section, &numFields);PYLITH_CHECK_ERROR(err);

  for (int iPoint=0; iPoint < numPoints; ++iPoint) {
    const PetscInt point = _points[iPoint];

    // Get list of currently constrained DOF
    PetscInt cdof;
    const PetscInt *cInd;
    err = PetscSectionGetConstraintDof(section, point, &cdof);PYLITH_CHECK_ERROR(err);
    err = PetscSectionGetConstraintIndices(section, point, &cInd);PYLITH_CHECK_ERROR(err);

    // Create array holding all constrained DOF
    int_array allCInd(cInd, cdof);

    // Verify other BC has not already constrained DOF
    const int numPrevious = _offsetLocal[iPoint];
    for (int iDOF=0; iDOF < numPrevious; ++iDOF)
      for (int jDOF=0; jDOF < numFixedDOF; ++jDOF)
        if (allCInd[iDOF] == _bcDOF[jDOF]) {
          std::ostringstream msg;
          msg << "Found multiple constraints on degrees of freedom at "
              << "point while setting up constraints for DirichletBC "
              << "boundary condition '" << _label << "'. "
	      << "Degree of freedom " << _bcDOF[jDOF] 
              << " is already constrained by another Dirichlet BC.";
          throw std::runtime_error(msg.str());
        } // if

    // Add in the ones for this DirichletBC BC
    for (int iDOF=0; iDOF < numFixedDOF; ++iDOF) {
      assert(_offsetLocal[iPoint]+iDOF < cdof);
      allCInd[_offsetLocal[iPoint]+iDOF] = _bcDOF[iDOF];
    } // for

    // Fill in rest of values not yet set (will be set by another DirichletBC BC)
    for (int iDOF=_offsetLocal[iPoint]+numFixedDOF; iDOF < cdof; ++iDOF) {
      allCInd[iDOF] = 999;
    } // for

    // Sort list of constrained DOF
    //   I need these sorted for my update algorithms to work properly
    std::sort(&allCInd[0], &allCInd[cdof]);

    // Update list of constrained DOF
    err = PetscSectionSetConstraintIndices(section, point, &allCInd[0]);PYLITH_CHECK_ERROR(err);
    if (numFields) {err = PetscSectionSetFieldConstraintIndices(section, point, 0, &allCInd[0]);PYLITH_CHECK_ERROR(err);}
  } // for

  PYLITH_METHOD_END;
} // setConstraints

// ----------------------------------------------------------------------
// Set values in field.
void
pylith::bc::DirichletBC::setField(const PylithScalar t,
				  const topology::Field& field)
{ // setField
  PYLITH_METHOD_BEGIN;

  const int numFixedDOF = _bcDOF.size();
  if (0 == numFixedDOF)
    PYLITH_METHOD_END;

  // Calculate spatial and temporal variation of value for BC.
  _calculateValue(t);

  // Get sections
  assert(_parameters);
  topology::Field& valueField = _parameters->get("value");
  topology::VecVisitorMesh valueVisitor(valueField);
  PetscScalar* valueArray = valueVisitor.localArray();

  topology::VecVisitorMesh fieldVisitor(field);
  PetscScalar* fieldArray = fieldVisitor.localArray();

  const int numPoints = _points.size();
  for (int iPoint=0; iPoint < numPoints; ++iPoint) {
    PetscInt p_bc = _points[iPoint];

    const PetscInt off = fieldVisitor.sectionOffset(p_bc);
    const PetscInt voff = valueVisitor.sectionOffset(p_bc);
    assert(numFixedDOF == valueVisitor.sectionDof(p_bc));

    for(PetscInt iDOF = 0; iDOF < numFixedDOF; ++iDOF) {
      assert(_bcDOF[iDOF] < fieldVisitor.sectionDof(p_bc));
      fieldArray[_bcDOF[iDOF]+off] = valueArray[voff+iDOF];
    } // for
  } // for

  PYLITH_METHOD_END;
} // setField

// ----------------------------------------------------------------------
// Set increment in values from t0 to t1 in field.
void
pylith::bc::DirichletBC::setFieldIncr(const PylithScalar t0,
				      const PylithScalar t1,
				      const topology::Field& field)
{ // setFieldIncr
  PYLITH_METHOD_BEGIN;

  const int numFixedDOF = _bcDOF.size();
  if (0 == numFixedDOF)
    PYLITH_METHOD_END;

  // Calculate spatial and temporal variation of value for BC.
  _calculateValueIncr(t0, t1);

  // Get sections
  assert(_parameters);
  topology::Field& valueField = _parameters->get("value");
  topology::VecVisitorMesh valueVisitor(valueField);
  PetscScalar* valueArray = valueVisitor.localArray();

  topology::VecVisitorMesh fieldVisitor(field);
  PetscScalar* fieldArray = fieldVisitor.localArray();

  const int numPoints = _points.size();
  for (int iPoint=0; iPoint < numPoints; ++iPoint) {
    PetscInt p_bc = _points[iPoint];

    const PetscInt off = fieldVisitor.sectionOffset(p_bc);
    const PetscInt voff = valueVisitor.sectionOffset(p_bc);
    assert(numFixedDOF == valueVisitor.sectionDof(p_bc));
    for(PetscInt iDOF=0; iDOF < numFixedDOF; ++iDOF) {
      assert(_bcDOF[iDOF] < fieldVisitor.sectionDof(p_bc));
      fieldArray[_bcDOF[iDOF]+off] = valueArray[voff+iDOF];
    } // for
  } // for

  PYLITH_METHOD_END;
} // setFieldIncr

// ----------------------------------------------------------------------
// Verify configuration is acceptable.
void
pylith::bc::DirichletBC::verifyConfiguration(const topology::Mesh& mesh) const
{ // verifyConfiguration
  PYLITH_METHOD_BEGIN;

  BoundaryCondition::verifyConfiguration(mesh);
  TimeDependent::verifyConfiguration(mesh);

  PYLITH_METHOD_END;
} // verifyConfiguration


// End of file 
