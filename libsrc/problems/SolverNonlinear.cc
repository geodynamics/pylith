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
  if (0 != _snes) {
    PetscErrorCode err = SNESDestroy(_snes); _snes = 0;
    CHECK_PETSC_ERROR(err);
  } // if
} // destructor

// ----------------------------------------------------------------------
// Initialize solver.
void
pylith::problems::SolverNonlinear::initialize(
			             topology::SolutionFields* fields,
				     const topology::Jacobian& jacobian,
				     Formulation::ResidualFn* residualFn,
				     Formulation::ArgsResidual* argsResidual,
				     Formulation::JacobianFn* jacobianFn,
				     Formulation::ArgsJacobian* argsJacobian)
{ // initialize
  assert(0 != fields);

  Solver::initialize(fields);

  PetscErrorCode err = 0;
  if (0 != _snes) {
    err = SNESDestroy(_snes); _snes = 0;
    CHECK_PETSC_ERROR(err);
  } // if    
  err = SNESCreate(fields->mesh().comm(), &_snes); CHECK_PETSC_ERROR(err);
  err = SNESSetFromOptions(_snes); CHECK_PETSC_ERROR(err);

  const topology::Field<topology::Mesh>& residual = fields->residual();
  const PetscVec residualVec = residual.vector();
  err = SNESSetFunction(_snes, residualVec, residualFn, argsResidual);
  CHECK_PETSC_ERROR(err);

  const PetscMat jacobianMat = jacobian.matrix();
  err = SNESSetJacobian(_snes, jacobianMat, jacobianMat,
			jacobianFn, argsJacobian);
  CHECK_PETSC_ERROR(err);
} // initialize

// ----------------------------------------------------------------------
// Solve the system.
void
pylith::problems::SolverNonlinear::solve(
			      topology::Field<topology::Mesh>* solution,
			      const topology::Jacobian& jacobian,
			      const topology::Field<topology::Mesh>& residual)
{ // solve
  assert(0 != solution);

  PetscErrorCode err = 0;

  // Update PetscVector view of field.
  residual.scatterSectionToVector();

  const PetscVec solutionVec = solution->vector();
  err = SNESSolve(_ksp, residualVec, solutionVec); CHECK_PETSC_ERROR(err);

  // Update section view of field.
  solution->scatterVectorToSection();
} // solve


// End of file
