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

#include "DirichletBC.hh" // implementation of object methods

#include "pylith/topology/Mesh.hh" // USES Mesh
#include "pylith/topology/Field.hh" // USES Field
#include "pylith/topology/Fields.hh" // USES Fields
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
} // destructor

// ----------------------------------------------------------------------
// Deallocate PETSc and local data structures.
void
pylith::bc::DirichletBC::deallocate(void)
{ // deallocate
  TimeDependentPoints::deallocate();
  feassemble::Constraint::deallocate();
} // deallocate
  
// ----------------------------------------------------------------------
// Initialize boundary condition.
void
pylith::bc::DirichletBC::initialize(const topology::Mesh& mesh,
				    const PylithScalar upDir[3])
{ // initialize
  if (0 == _bcDOF.size())
    return;

  _getPoints(mesh);

  assert(_normalizer);
  const PylithScalar lengthScale = _normalizer->lengthScale();
  _queryDatabases(mesh, lengthScale, "displacement");
} // initialize

// ----------------------------------------------------------------------
// Set number of degrees of freedom that are constrained at points in field.
void
pylith::bc::DirichletBC::setConstraintSizes(const topology::Field<topology::Mesh>& field)
{ // setConstraintSizes
  const int numFixedDOF = _bcDOF.size();
  if (0 == numFixedDOF)
    return;

  PetscSection section = field.petscSection();assert(section);
  PetscInt numFields;
  PetscErrorCode err = 0;

  // Set constraints in field
  const int numPoints = _points.size();
  _offsetLocal.resize(numPoints);
  err = PetscSectionGetNumFields(section, &numFields);CHECK_PETSC_ERROR(err);

  for (int iPoint=0; iPoint < numPoints; ++iPoint) {
    const PetscInt point = _points[iPoint];

    PetscInt dof, cdof;
    err = PetscSectionGetDof(section, point, &dof);CHECK_PETSC_ERROR(err);
    err = PetscSectionGetConstraintDof(section, point, &cdof);CHECK_PETSC_ERROR(err);
    if (cdof + numFixedDOF > dof) {
      std::ostringstream msg;
      msg
	<< "Found overly constrained point while setting up constraints for\n"
	<< "DirichletBC boundary condition '" << _label << "'.\n"
	<< "Number of DOF at point " << point << " is " << dof
	<< "\nand number of attempted constraints is " << cdof+numFixedDOF << ".";
      throw std::runtime_error(msg.str());
    } // if
    _offsetLocal[iPoint] = cdof;
    err = PetscSectionAddConstraintDof(section, point, numFixedDOF);CHECK_PETSC_ERROR(err);
    // We should be specifying what field the BC is for
    if (numFields) {err = PetscSectionAddFieldConstraintDof(section, point, 0, numFixedDOF);CHECK_PETSC_ERROR(err);}
  } // for
} // setConstraintSizes

// ----------------------------------------------------------------------
// Set which degrees of freedom are constrained at points in field.
void
pylith::bc::DirichletBC::setConstraints(const topology::Field<topology::Mesh>& field)
{ // setConstraints
  const int numFixedDOF = _bcDOF.size();
  if (0 == numFixedDOF)
    return;

  PetscSection section = field.petscSection();assert(section);
  PetscInt numFields;
  PetscErrorCode err = 0;

  const int numPoints = _points.size();
  err = PetscSectionGetNumFields(section, &numFields);CHECK_PETSC_ERROR(err);

  for (int iPoint=0; iPoint < numPoints; ++iPoint) {
    const PetscInt point = _points[iPoint];

    // Get list of currently constrained DOF
    PetscInt cdof;
    const PetscInt *cInd;
    err = PetscSectionGetConstraintDof(section, point, &cdof);CHECK_PETSC_ERROR(err);
    err = PetscSectionGetConstraintIndices(section, point, &cInd);CHECK_PETSC_ERROR(err);

    // Create array holding all constrained DOF
    int_array allCInd(cInd, cdof);

    // Verify other BC has not already constrained DOF
    const int numPrevious = _offsetLocal[iPoint];
    for (int iDOF=0; iDOF < numPrevious; ++iDOF)
      for (int jDOF=0; jDOF < numFixedDOF; ++jDOF)
        if (allCInd[iDOF] == _bcDOF[jDOF]) {
          std::ostringstream msg;
          msg << "Found multiple constraints on degrees of freedom at\n"
              << "point while setting up constraints for DirichletBC\n"
              << "boundary condition '" << _label << "'.\n"
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
    err = PetscSectionSetConstraintIndices(section, point, &allCInd[0]);CHECK_PETSC_ERROR(err);
    if (numFields) {err = PetscSectionSetFieldConstraintIndices(section, point, 0, &allCInd[0]);CHECK_PETSC_ERROR(err);}
  } // for
} // setConstraints

// ----------------------------------------------------------------------
// Set values in field.
void
pylith::bc::DirichletBC::setField(const PylithScalar t,
				  const topology::Field<topology::Mesh>& field)
{ // setField
  const int numFixedDOF = _bcDOF.size();
  if (0 == numFixedDOF)
    return;

  // Calculate spatial and temporal variation of value for BC.
  _calculateValue(t);

  assert(_parameters);
  topology::Field<topology::Mesh>& valueField = _parameters->get("value");
  PetscSection fieldSection = field.petscSection();assert(fieldSection);

  PetscScalar* valueArray = valueField.getLocalArray();
  PetscScalar* fieldArray = field.getLocalArray();

  const int numPoints = _points.size();
  for (int iPoint=0; iPoint < numPoints; ++iPoint) {
    PetscInt p_bc = _points[iPoint];

    const PetscInt off = field.sectionOffset(p_bc);
    const PetscInt voff = valueField.sectionOffset(p_bc);
    assert(numFixedDOF == valueField.sectionDof(p_bc));

    for(PetscInt iDOF = 0; iDOF < numFixedDOF; ++iDOF) {
      assert(_bcDOF[iDOF] < field.sectionDof(p_bc));
      fieldArray[_bcDOF[iDOF]+off] = valueArray[voff+iDOF];
    } // for
  } // for
  valueField.restoreLocalArray(&valueArray);
  field.restoreLocalArray(&fieldArray);
} // setField

// ----------------------------------------------------------------------
// Set increment in values from t0 to t1 in field.
void
pylith::bc::DirichletBC::setFieldIncr(const PylithScalar t0,
				      const PylithScalar t1,
				      const topology::Field<topology::Mesh>& field)
{ // setFieldIncr
  const int numFixedDOF = _bcDOF.size();
  if (0 == numFixedDOF)
    return;

  // Calculate spatial and temporal variation of value for BC.
  _calculateValueIncr(t0, t1);

  assert(_parameters);
  topology::Field<topology::Mesh>& valueField = _parameters->get("value");
  PetscSection fieldSection = field.petscSection();assert(fieldSection);

  PetscScalar* valueArray = valueField.getLocalArray();
  PetscScalar* fieldArray = field.getLocalArray();

  const int numPoints = _points.size();
  for (int iPoint=0; iPoint < numPoints; ++iPoint) {
    PetscInt p_bc = _points[iPoint];

    const PetscInt off = field.sectionOffset(p_bc);
    const PetscInt voff = valueField.sectionOffset(p_bc);
    assert(numFixedDOF == valueField.sectionDof(p_bc));
    for(PetscInt iDOF=0; iDOF < numFixedDOF; ++iDOF) {
      assert(_bcDOF[iDOF] < field.sectionDof(p_bc));
      fieldArray[_bcDOF[iDOF]+off] = valueArray[voff+iDOF];
    } // for
  } // for
  valueField.restoreLocalArray(&valueArray);
  field.restoreLocalArray(&fieldArray);
} // setFieldIncr

// ----------------------------------------------------------------------
// Verify configuration is acceptable.
void
pylith::bc::DirichletBC::verifyConfiguration(const topology::Mesh& mesh) const
{ // verifyConfiguration
  BoundaryCondition::verifyConfiguration(mesh);
  TimeDependent::verifyConfiguration(mesh);
} // verifyConfiguration


// End of file 
