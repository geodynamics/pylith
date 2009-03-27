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
#include "pylith/topology/Jacobian.hh" // USES Jacobian

#include <petscksp.h> // USES PetscKSP

#include "pylith/utils/petscerror.h" // USES CHECK_PETSC_ERROR

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
  deallocate();
} // destructor

// ----------------------------------------------------------------------
// Deallocate data structures.
void
pylith::problems::SolverLinear::deallocate(void)
{ // deallocate
  if (0 != _ksp) {
    PetscErrorCode err = KSPDestroy(_ksp); _ksp = 0;
    CHECK_PETSC_ERROR(err);
    std::cout << "DESTROYED KSP." << std::endl;
  } // if
} // deallocate
  
// ----------------------------------------------------------------------
// Set initial guess zero flag.
void
pylith::problems::SolverLinear::initialGuessZero(const bool value)
{ // initialGuessZero
  assert(0 != _ksp);

  PetscTruth flag = (value) ? PETSC_FALSE : PETSC_TRUE;
  PetscErrorCode err = KSPSetInitialGuessNonzero(_ksp, flag);
  CHECK_PETSC_ERROR(err);
} // initialGuessZero

// ----------------------------------------------------------------------
// Initialize solver.
void
pylith::problems::SolverLinear::initialize(
				   const topology::SolutionFields& fields,
				   const topology::Jacobian& jacobian,
				   Formulation* formulation)
{ // initialize
  assert(0 != formulation);

  Solver::initialize(fields, jacobian, formulation);

  PetscErrorCode err = 0;
  if (0 != _ksp) {
    err = KSPDestroy(_ksp); _ksp = 0;
    CHECK_PETSC_ERROR(err);
  } // if    
  err = KSPCreate(fields.mesh().comm(), &_ksp); CHECK_PETSC_ERROR(err);
  err = KSPSetFromOptions(_ksp); CHECK_PETSC_ERROR(err);
} // initialize

// ----------------------------------------------------------------------
// Solve the system.
void
pylith::problems::SolverLinear::solve(
			      topology::Field<topology::Mesh>* solution,
			      const topology::Jacobian& jacobian,
			      const topology::Field<topology::Mesh>& residual)
{ // solve
  assert(0 != solution);

  PetscErrorCode err = 0;

  // Update PetscVector view of field.
  residual.scatterSectionToVector();

  const PetscMat jacobianMat = jacobian.matrix();
  err = KSPSetOperators(_ksp, jacobianMat, jacobianMat, 
			DIFFERENT_NONZERO_PATTERN); CHECK_PETSC_ERROR(err);

  const PetscVec residualVec = residual.vector();
  const PetscVec solutionVec = solution->vector();
  err = KSPSolve(_ksp, residualVec, solutionVec); CHECK_PETSC_ERROR(err);

  // Update section view of field.
  solution->scatterVectorToSection();
} // solve


// End of file
