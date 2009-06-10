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

//#define FIELD_SPLIT

#if defined(FIELD_SPLIT)
#include <petscmesh_solvers.hh> // USES constructFieldSplit()
#endif

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
  } // if
} // deallocate
  
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
  err = KSPSetInitialGuessNonzero(_ksp, PETSC_FALSE); CHECK_PETSC_ERROR(err);
  err = KSPSetFromOptions(_ksp); CHECK_PETSC_ERROR(err);

  const topology::Field<topology::Mesh>& residual = fields.get("residual");

  // Check for fibration
  if (residual.section()->getNumSpaces() > 0) {
    PC pc;

    err = KSPGetPC(_ksp, &pc); CHECK_PETSC_ERROR(err);
    err = PCSetType(pc, PCFIELDSPLIT); CHECK_PETSC_ERROR(err);

#if defined(FIELD_SPLIT)
    constructFieldSplit(residual.section(), pc);
#endif
  }
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
