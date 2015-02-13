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
// Copyright (c) 2010-2015 University of California, Davis
//
// See COPYING for license information.
//
// ======================================================================
//

#include <portinfo>

#include "SolverLinear.hh" // implementation of class methods

#include "Formulation.hh" // USES Formulation

#include "pylith/topology/Mesh.hh" // USES Mesh
#include "pylith/topology/Field.hh" // USES Field
#include "pylith/topology/SolutionFields.hh" // USES SolutionFields
#include "pylith/topology/Jacobian.hh" // USES Jacobian

#include "pylith/utils/EventLogger.hh" // USES EventLogger

#include <petscksp.h> // USES PetscKSP

#include "pylith/utils/error.h" // USES PYLITH_CHECK_ERROR

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

  PetscErrorCode err = KSPDestroy(&_ksp);PYLITH_CHECK_ERROR(err);

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
  err = KSPDestroy(&_ksp);PYLITH_CHECK_ERROR(err);
  err = KSPCreate(fields.mesh().comm(), &_ksp);PYLITH_CHECK_ERROR(err);
  err = KSPSetInitialGuessNonzero(_ksp, PETSC_FALSE);PYLITH_CHECK_ERROR(err);
  err = KSPSetFromOptions(_ksp);PYLITH_CHECK_ERROR(err);

  if (formulation->splitFields()) {
    PetscPC pc = 0;
    err = KSPGetPC(_ksp, &pc);PYLITH_CHECK_ERROR(err);
    _setupFieldSplit(&pc, formulation, jacobian, fields);
  } // if

  _createNullSpace(fields);

  PYLITH_METHOD_END;
} // initialize

// ----------------------------------------------------------------------
// Solve the system.
void
pylith::problems::SolverLinear::solve(topology::Field* solution,
				      topology::Jacobian* jacobian,
				      const topology::Field& residual)
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
  residual.scatterLocalToGlobal();

  _logger->eventEnd(scatterEvent);
  _logger->eventBegin(setupEvent);

  PetscErrorCode err = 0;
  const PetscMat jacobianMat = jacobian->matrix();
  err = KSPSetOperators(_ksp, jacobianMat, jacobianMat);PYLITH_CHECK_ERROR(err);
  jacobian->resetValuesChanged();

  const PetscVec residualVec = residual.globalVector();
  const PetscVec solutionVec = solution->globalVector();

  _logger->eventEnd(setupEvent);
  _logger->eventBegin(solveEvent);

  err = KSPSolve(_ksp, residualVec, solutionVec); PYLITH_CHECK_ERROR(err);

  _logger->eventEnd(solveEvent);
  _logger->eventBegin(scatterEvent);

  // Update section view of field.
  solution->scatterGlobalToLocal();

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
