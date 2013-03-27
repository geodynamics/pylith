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
// Copyright (c) 2010-2012 University of California, Davis
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
  PYLITH_METHOD_BEGIN;

  Solver::deallocate();

  PetscErrorCode err = KSPDestroy(&_ksp);CHECK_PETSC_ERROR(err);

  PYLITH_METHOD_END;
} // deallocate
  
// ----------------------------------------------------------------------
// Initialize solver.
void
pylith::problems::SolverLinear::initialize(const topology::SolutionFields& fields,
					   const topology::Jacobian& jacobian,
					   Formulation* formulation)
{ // initialize
  PYLITH_METHOD_BEGIN;

  assert(formulation);

  _initializeLogger();
  Solver::initialize(fields, jacobian, formulation);

  PetscErrorCode err = 0;
  err = KSPDestroy(&_ksp);CHECK_PETSC_ERROR(err);
  err = KSPCreate(fields.mesh().comm(), &_ksp);CHECK_PETSC_ERROR(err);
  err = KSPSetInitialGuessNonzero(_ksp, PETSC_FALSE);CHECK_PETSC_ERROR(err);
  err = KSPSetFromOptions(_ksp);CHECK_PETSC_ERROR(err);

  if (formulation->splitFields()) {
    PetscPC pc = 0;
    err = KSPGetPC(_ksp, &pc);CHECK_PETSC_ERROR(err);
    _setupFieldSplit(&pc, formulation, jacobian, fields);
  } // if

  PYLITH_METHOD_END;
} // initialize

// ----------------------------------------------------------------------
// Solve the system.
void
pylith::problems::SolverLinear::solve(topology::Field<topology::Mesh>* solution,
				      topology::Jacobian* jacobian,
				      const topology::Field<topology::Mesh>& residual)
{ // solve
  PYLITH_METHOD_BEGIN;

  assert(solution);
  assert(jacobian);
  assert(_formulation);

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
    err = KSPSetOperators(_ksp, jacobianMat, jacobianMat, SAME_PRECONDITIONER);CHECK_PETSC_ERROR(err);
  } else {
    err = KSPSetOperators(_ksp, jacobianMat, jacobianMat, DIFFERENT_NONZERO_PATTERN);CHECK_PETSC_ERROR(err);
  } // else
  jacobian->resetValuesChanged();

  const PetscVec residualVec = residual.globalVector();
  const PetscVec solutionVec = solution->globalVector();

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

  PYLITH_METHOD_END;
} // solve

// ----------------------------------------------------------------------
// Initialize logger.
void
pylith::problems::SolverLinear::_initializeLogger(void)
{ // initializeLogger
  PYLITH_METHOD_BEGIN;

  delete _logger; _logger = new utils::EventLogger;assert(_logger);
  _logger->className("SolverLinear");
  _logger->initialize();
  _logger->registerEvent("SoLi setup");
  _logger->registerEvent("SoLi solve");
  _logger->registerEvent("SoLi scatter");

  PYLITH_METHOD_END;
} // initializeLogger


// End of file
