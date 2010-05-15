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
#include "pylith/problems/Formulation.hh" // USES Formulation
#include "pylith/utils/EventLogger.hh" // USES EventLogger

#include <petscksp.h> // USES PetscKSP

#include "pylith/utils/petscerror.h" // USES CHECK_PETSC_ERROR

#define FIELD_SPLIT
#define NEW_FAULT_PRECONDITIONER

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
  Solver::deallocate();

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

  _initializeLogger();
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
    const ALE::Obj<topology::Mesh::SieveMesh>& sieveMesh = fields.mesh().sieveMesh();
    PC pc;

    err = KSPGetPC(_ksp, &pc); CHECK_PETSC_ERROR(err);
    err = PCSetType(pc, PCFIELDSPLIT); CHECK_PETSC_ERROR(err);
    err = PCSetOptionsPrefix(pc, "fs_"); CHECK_PETSC_ERROR(err);
    err = PCSetFromOptions(pc); CHECK_PETSC_ERROR(err);
#if defined(FIELD_SPLIT)
    constructFieldSplit(residual.section(), sieveMesh->getFactory()->getGlobalOrder(sieveMesh, "default", residual.section()), residual.vector(), pc);
#if defined(NEW_FAULT_PRECONDITIONER)
    if (residual.section()->getNumSpaces() > sieveMesh->getDimension()) {
      KSP     *ksps;
      Mat      A, M;
      PetscInt num, m, n;

      err = PCFieldSplitGetSubKSP(pc, &num, &ksps); CHECK_PETSC_ERROR(err);
      // Put in PC matrix for fault
      MatStructure flag;
      err = KSPGetOperators(ksps[num-1], &A, PETSC_NULL, &flag); CHECK_PETSC_ERROR(err);
      err = PetscObjectReference((PetscObject) A); CHECK_PETSC_ERROR(err);
      err = MatGetLocalSize(A, &m, &n); CHECK_PETSC_ERROR(err);
      err = MatCreate(sieveMesh->comm(), &M); CHECK_PETSC_ERROR(err);
      err = MatSetSizes(M, m, n, PETSC_DECIDE, PETSC_DECIDE); CHECK_PETSC_ERROR(err);
      err = MatSeqAIJSetPreallocation(M, 1, PETSC_NULL); CHECK_PETSC_ERROR(err);
      err = MatMPIAIJSetPreallocation(M, 1, PETSC_NULL, 0, PETSC_NULL); CHECK_PETSC_ERROR(err);
      err = MatSetFromOptions(M); CHECK_PETSC_ERROR(err);
      err = KSPSetOperators(ksps[num-1], A, M, flag); CHECK_PETSC_ERROR(err);
      // Create a mapping to indices for that space (might be in FS)
      err = PetscFree(ksps); CHECK_PETSC_ERROR(err);
    } // if
#endif
#endif
  }
} // initialize

// ----------------------------------------------------------------------
// Solve the system.
void
pylith::problems::SolverLinear::solve(
			      topology::Field<topology::Mesh>* solution,
			      topology::Jacobian* jacobian,
			      const topology::Field<topology::Mesh>& residual)
{ // solve
  assert(0 != solution);
  assert(0 != jacobian);
  assert(0 != _formulation);

  const int setupEvent = _logger->eventId("SoLi setup");
  const int solveEvent = _logger->eventId("SoLi solve");
  const int scatterEvent = _logger->eventId("SoLi scatter");
  _logger->eventBegin(scatterEvent);

  // Update PetscVector view of field.
  residual.scatterSectionToVector();

  _logger->eventEnd(scatterEvent);
  _logger->eventBegin(setupEvent);

  PetscErrorCode err = 0;
  const PetscMat jacobianMat = jacobian->matrix();
  if (!jacobian->valuesChanged()) {
    err = KSPSetOperators(_ksp, jacobianMat, jacobianMat, 
			  SAME_PRECONDITIONER); CHECK_PETSC_ERROR(err);
  } else {
    err = KSPSetOperators(_ksp, jacobianMat, jacobianMat, 
			  DIFFERENT_NONZERO_PATTERN); CHECK_PETSC_ERROR(err);
  } // else
  jacobian->resetValuesChanged();

  const PetscVec residualVec = residual.vector();
  const PetscVec solutionVec = solution->vector();

  _logger->eventEnd(setupEvent);
  _logger->eventBegin(solveEvent);

  err = KSPSolve(_ksp, residualVec, solutionVec); CHECK_PETSC_ERROR(err);

  _logger->eventEnd(solveEvent);
  _logger->eventBegin(scatterEvent);

  // Update section view of field.
  solution->scatterVectorToSection();

  _logger->eventEnd(scatterEvent);

  // Update rate fields to be consistent with current solution.
  _formulation->calcRateFields();
} // solve

// ----------------------------------------------------------------------
// Initialize logger.
void
pylith::problems::SolverLinear::_initializeLogger(void)
{ // initializeLogger
  delete _logger; _logger = new utils::EventLogger;
  assert(0 != _logger);
  _logger->className("SolverLinear");
  _logger->initialize();
  _logger->registerEvent("SoLi setup");
  _logger->registerEvent("SoLi solve");
  _logger->registerEvent("SoLi scatter");
} // initializeLogger


// End of file
