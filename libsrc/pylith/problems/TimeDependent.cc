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

#include "TimeDependent.hh" // implementation of class methods

// ----------------------------------------------------------------------
// Constructor
pylith::problems::TimeDependent::TimeDependent(void)
{ // constructor
} // constructor

// ----------------------------------------------------------------------
// Destructor
pylith::problems::TimeDependent::~TimeDependent(void)
{ // destructor
  deallocate();
} // destructor

// ----------------------------------------------------------------------
// Deallocate PETSc and local data structures.
void
pylith::problems::TimeDependent::deallocate(void)
{ // deallocate
  PYLITH_METHOD_BEGIN;

  PetscErrorCode err = TSDestroy(&_ts);PYLITH_CHECK_ERROR(_ts);

  PYLITH_METHOD_END;
} // deallocate

// ----------------------------------------------------------------------
// Set problem type.
void
pylith::problems::TimeDependent::problemType(const ProblemTypeEnum value)
{ // problemType
  _problemType = value;
} // problemType

// ----------------------------------------------------------------------
// Get problem type.
pylith::problems::TimeDependent::ProblemTypeEnum
pylith::problems::TimeDependent::problemType(void) const
{ // problemType
  return _problemType;
} // problemType

// ----------------------------------------------------------------------
// Set start time for problem.
void
pylith::problems::TimeDependent::startTime(const PetscReal value)
{ // startTime
  _startTime = value;
} // startTime

// ----------------------------------------------------------------------
// Get start time for problem.
PetscReal
pylith::problems::TimeDependent::startTime(void) const
{ // startTime
  return _startTime;
} // startTime

// ----------------------------------------------------------------------
// Set total time for problem.
void
pylith::problems::TimeDependent::totalTime(const PetscReal value)
{ // totalTime
  PYLITH_METHOD_BEGIN;

  if (value < 0.0) {
    std::ostringstream msg;
    msg << "Total time (nondimensional) for problem (" << value << ") must be positive.";
    throw std::runtime_error(msg.str());
  } // if
  _totalTime = value;

  PYLITH_METHOD_END;
} // totalTime

// ----------------------------------------------------------------------
// Get total time for problem.
PetscReal
pylith::problems::TimeDependent::totalTime(void) const
{ // totalTime
  return _totalTime;
} // totalTime

// ----------------------------------------------------------------------
// Set initial time step for problem.
void
pylith::problems::TimeDependent::dtInitial(const PetscReal value)
{ // dtInitial
  PYLITH_METHOD_BEGIN;

  if (value < 0.0) {
    std::ostringstream msg;
    msg << "Initial time step (nondimensional) for problem (" << value << ") must be positive.";
    throw std::runtime_error(msg.str());
  } // if
  _dtInitial = value;

  PYLITH_METHOD_END;
} // dtInitial

// ----------------------------------------------------------------------
// Get initial time step for problem.
PetscReal
pylith::problems::TimeDependent::dtInitial(void) const
{ // dtInitial
  return _dtInitial;
} // dtInitial

// ----------------------------------------------------------------------
// Initialize.
void
pylith::problems::TimeDependent::initialize(pylith::topology::Field* solution,
					    pylith::topology::Jacobian* jacobianRHS)
{ // initialize
  PYLITH_METHOD_BEGIN;

  assert(solution);
  assert(jacobianRHS);

  _solution = solution;
  _jacobianRHS = jacobianRHS;

  PetscErrorCode err = TSDestroy(&_ts);PYLITH_CHECK_ERROR(err);assert(!_ts);
  err = TSCreate(comm, &_ts);PYLITH_CHECK_ERROR(err);assert(_ts);
  
  // Set time stepping paramters.
  switch (_problemType) {
  case LINEAR:
    err = TSSetProblemType(_ts, TS_LINEAR);PYLITH_CHECK_ERROR(err);
    break;
  case NONLINEAR:
    err = TSSetProblemType(_ts, TS_NONLINEAR);PYLITH_CHECK_ERROR(err);
    break;
  default:
    assert(0);
    throw std::logic_error("Unknown problem type.");
  } // switch
  err = TSSetInitialTimeStep(_ts, _initialDt, _startTime);PYLITH_CHECK_ERROR(err);
  err = TSSetDuration(_ts, _totalTIme);PYLITH_CHECK_ERROR(err);

  // Set solution.
  err = TSSetSolution(_ts, _solution->globalVector());PYLITH_CHECKE_ERROR(err);

  // Set callbacks.
  err = TSSetPreStep(_ts, prestep);PYLITH_CHECK_ERROR(err);
  err = TSSetPreStep(_ts, poststep);PYLITH_CHECK_ERROR(err);
  err = TSSetRHSJacoabian(_ts, _jacobianRHS->matrix(), precondMatrix, reformRHSJacobian, (void*)this);PYLITH_CHECK_ERROR(err);
  err = TSSetRHSFunction(_ts, residual, reformResidual, (void*)this);PYLITH_CHECK_ERROR(err);

  // :TODO: Implement setting initial conditions.

  // Setup time stepper.
  err = TSSetUp(_ts);PYLITH_CHECK_ERROR(err);

  PYLITH_METHOD_END;
} // initialize

// ----------------------------------------------------------------------
// Solve time dependent problem.
void
pylith::problems::TimeDependent::solve(void)
{ // solve
  PYLITH_METHOD_BEGIN;

  PetscErrorCode err = TSSolve(_ts, NULL);PYLITH_CHECK_ERROR(err);

  PYLITH_METHOD_END;
} // solve


// ----------------------------------------------------------------------
// Callback method for operations before advancing solution one time step.
PetscErrorCode
pylith::problems::TimeDependent::prestep(PetscTS ts)
{ // prestep
  PYLITH_METHOD_BEGIN;

  TimeDependent* problem = NULL;
  PetscErrorCode err = TSGetApplicationContext(ts, (void*)problem);PYLITH_CHECK_ERROR(err);assert(problem);
  problem->prestep();

  PYLITH_METHOD_RETURN(0);
} // prestep


// ----------------------------------------------------------------------
// Callback method for operations after advancing solution one time step.
void
pylith::problems::TimeDependent::poststep(PetscTS ts)
{ // poststep
  PYLITH_METHOD_BEGIN;

  TimeDependent* problem = NULL;
  PetscErrorCode err = TSGetApplicationContext(ts, (void*)problem);PYLITH_CHECK_ERROR(err);assert(problem);
  problem->poststep();

  PYLITH_METHOD_RETURN(0);
} // poststep


// ----------------------------------------------------------------------
// Perform operations before advancing solution one time step.
void
pylith::problems::TimeDependent::prestep(void)
{ // prestep
  PYLITH_METHOD_BEGIN;

  // Set constraints.

  

  PYLITH_METHOD_END;
} // prestep

// ----------------------------------------------------------------------
// Perform operations after advancing solution one time step.
void
pylith::problems::TimeDependent::poststep(void)
{ // poststep
  PYLITH_METHOD_BEGIN;

  PYLITH_METHOD_END;
} // poststep

// End of file
