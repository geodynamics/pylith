// -*- C++ -*-
//
// ======================================================================
//
//                           Brad T. Aagaard
//                        U.S. Geological Survey
//
// {LicenseText}
//
// ======================================================================
//

#include <portinfo>

#include "SolverNonlinear.hh" // implementation of class methods

#include "Formulation.hh" // USES Formulation
#include "pylith/topology/SolutionFields.hh" // USES SolutionFields
#include "pylith/topology/Jacobian.hh" // USES Jacobian

#include <petscsnes.h> // USES PetscSNES

#include "pylith/utils/petscerror.h" // USES CHECK_PETSC_ERROR

// ----------------------------------------------------------------------
// Constructor
pylith::problems::SolverNonlinear::SolverNonlinear(void) :
  _snes(0)
{ // constructor
} // constructor

// ----------------------------------------------------------------------
// Destructor
pylith::problems::SolverNonlinear::~SolverNonlinear(void)
{ // destructor
  deallocate();
} // destructor

// ----------------------------------------------------------------------
// Deallocate data structures.
void
pylith::problems::SolverNonlinear::deallocate(void)
{ // deallocate
  if (0 != _snes) {
    PetscErrorCode err = SNESDestroy(_snes); _snes = 0;
    CHECK_PETSC_ERROR(err);
  } // if
} // deallocate
  
// ----------------------------------------------------------------------
// Initialize solver.
void
pylith::problems::SolverNonlinear::initialize(
			             const topology::SolutionFields& fields,
				     const topology::Jacobian& jacobian,
				     Formulation* formulation)
{ // initialize
  assert(0 != formulation);

  Solver::initialize(fields, jacobian, formulation);

  PetscErrorCode err = 0;
  if (0 != _snes) {
    err = SNESDestroy(_snes); _snes = 0;
    CHECK_PETSC_ERROR(err);
  } // if    
  err = SNESCreate(fields.mesh().comm(), &_snes); CHECK_PETSC_ERROR(err);
  err = SNESSetFromOptions(_snes); CHECK_PETSC_ERROR(err);

  const topology::Field<topology::Mesh>& residual = fields.get("residual");
  const PetscVec residualVec = residual.vector();
  err = SNESSetFunction(_snes, residualVec, reformResidual,
			(void*) formulation);
  CHECK_PETSC_ERROR(err);

  PetscMat jacobianMat = jacobian.matrix();
  err = SNESSetJacobian(_snes, jacobianMat, jacobianMat, reformJacobian,
			(void*) formulation);
  CHECK_PETSC_ERROR(err);
} // initialize

// ----------------------------------------------------------------------
// Solve the system.
void
pylith::problems::SolverNonlinear::solve(
			      topology::Field<topology::Mesh>* solveSoln,
			      const topology::Jacobian& jacobian,
			      const topology::Field<topology::Mesh>& residual)
{ // solve
  assert(0 != solveSoln);

  PetscErrorCode err = 0;

  const PetscVec solveSolnVec = solveSoln->vector();
  err = SNESSolve(_snes, PETSC_NULL, solveSolnVec); CHECK_PETSC_ERROR(err);
  
  // Update section view of field.
  solveSoln->scatterVectorToSection();
} // solve

// ----------------------------------------------------------------------
// Generic C interface for reformResidual for integration with
// PETSc SNES solvers.
PetscErrorCode
pylith::problems::SolverNonlinear::reformResidual(PetscSNES snes,
						  PetscVec tmpSolutionVec,
						  PetscVec tmpResidualVec,
						  void* context)
{ // reformResidual
  assert(0 != context);
  Formulation* formulation = (Formulation*) context;
  assert(0 != formulation);

  // Reform residual
  formulation->reformResidual(&tmpResidualVec, &tmpSolutionVec);

  return 0;
} // reformResidual

// ----------------------------------------------------------------------
// Generic C interface for reformJacobian for integration with
// PETSc SNES solvers.
PetscErrorCode
pylith::problems::SolverNonlinear::reformJacobian(PetscSNES snes,
						  PetscVec tmpSolutionVec,
						  PetscMat* jacobianMat,
						  PetscMat* preconditionerMat,
						  MatStructure* preconditionerLayout,
						  void* context)
{ // reformJacobian
  assert(0 != context);
  Formulation* formulation = (Formulation*) context;
  assert(0 != formulation);

  formulation->reformJacobian(&tmpSolutionVec);

  return 0;
} // reformJacobian


// End of file
