// -*- C++ -*-
//
// ======================================================================
//
// Brad T. Aagaard, U.S. Geological Survey
// Charles A. Williams, GNS Science
// Matthew G. Knepley, University of Chicago
//
// This code was developed as part of the Computational Infrastructure
// for Geodynamics (http://geodynamics.org).
//
// Copyright (c) 2010 University of California, Davis
//
// See COPYING for license information.
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

// ----------------------------------------------------------------------
typedef pylith::topology::Mesh::SieveMesh SieveMesh;
typedef pylith::topology::Mesh::RealSection RealSection;

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

  if (formulation->splitFields()) {
    PetscPC pc = 0;
    err = KSPGetPC(_ksp, &pc); CHECK_PETSC_ERROR(err);
    _ctx.pc = pc;
    _setupFieldSplit(&pc, &_precondMatrix, formulation, fields);
  } // if
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

  // Update KSP operators with custom preconditioner if necessary.
  const ALE::Obj<RealSection>& solutionSection = solution->section();
  assert(!solutionSection.isNull());
  const ALE::Obj<SieveMesh>& sieveMesh = solution->mesh().sieveMesh();
  assert(!sieveMesh.isNull());
#if 0
  if (solutionSection->getNumSpaces() > sieveMesh->getDimension() &&
      0 != _precondMatrix) {
    
    PetscPC pc = 0;
    PetscKSP *ksps = 0;
    PetscMat A = 0;
    PetscInt num = 0;
    
    err = KSPSetUp(_ksp); CHECK_PETSC_ERROR(err);
    err = KSPGetPC(_ksp, &pc); CHECK_PETSC_ERROR(err);
    err = PCFieldSplitGetSubKSP(pc, &num, &ksps); CHECK_PETSC_ERROR(err);
    assert(solutionSection->getNumSpaces() == num);

    MatStructure flag;
    err = KSPGetOperators(ksps[num-1], &A, 
			  PETSC_NULL, &flag); CHECK_PETSC_ERROR(err);
    err = PetscObjectReference((PetscObject) A); CHECK_PETSC_ERROR(err);
    err = KSPSetOperators(ksps[num-1], A, _precondMatrix, 
			  flag); CHECK_PETSC_ERROR(err);
    err = PetscFree(ksps); CHECK_PETSC_ERROR(err);
  } // if
#endif

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
