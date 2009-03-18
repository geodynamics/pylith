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

#include "SolverLinear.hh" // implementation of class methods

#include "pylith/topology/SolutionFields.hh" // USES SolutionFields

#include <petscksp.h> // USES PetscKSP

// ----------------------------------------------------------------------
// Constructor
pylith::problems::SolverLinear::SolverLinear(void) :
  _ksp(0)
{ // constructor
} // constructor

// ----------------------------------------------------------------------
// Destructor
pylith::problems::SolverLinear::~SolverLinear(void)
{ // destructor
  if (0 != _ksp) {
    PetscErrorCode err = KSPDestroy(_ksp); _ksp = 0;
    CHECK_PETSC_ERROR(err);
  } // if
} // destructor

// ----------------------------------------------------------------------
// Set initial guess nonzero flag.
void
pylith::problems::SolverLinear::initialGuessNonzero(const bool value)
{ // initialGuessNonzero
  assert(0 != _ksp);

  PetscTruth flag = (value) ? PETSC_TRUE : PETSC_FALSE;
  KSPSetInitialGuessNonzero(_ksp, flag);
} // initialGuessNonzero

// ----------------------------------------------------------------------
// Initialize solver.
void
pylith::problems::SolverLinear::initialize(topology::SolutionFields* fields)
{ // initialize
  assert(0 != fields);

  Solver::initialize(fields);

  PetscErrorCode err = 0;
  if (0 != _ksp) {
    err = KSPDestroy(_ksp); _ksp = 0;
    CHECK_PETSC_ERROR(err);
  } // if    
  err = KSPCreate(fields->mesh().comm(), &_ksp); CHECK_PETSC_ERROR(err);
  err = KSPSetFromOptions(_ksp); CHECK_PETSC_ERROR(err);
} // initialize

// ----------------------------------------------------------------------
// Solve the system.
void
pylith::problems::SolverLinear::solve(topology::Field<topology::Mesh>* solution,
				      const topology::Jacobian& jacobian,
				      const topology::Field<topology::Mesh>& residual)
{ // solve
  assert(0 != solution);

  PetscErrorCode err = 0;

  // Update PetscVector view of field.
  residual->scatterSectionToVector();

  err = KSPSetOperators(_ksp, jacobianMat, jacobianMat, 
			DIFFERENT_NONZERO_PATTERN); CHECK_PETSC_ERROR(err);

  const PetscVec residualVec = residual.vector();
  const PetscVec solutionVec = solution->vector();
  err = KSPSolve(_ksp, residualVec, solutionVec); CHECK_PETSC_ERROR(err);

  // Update section view of field.
  solution->scatterVectorToSection();
} // solve


// End of file
