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
  err = TSSetFromOptions(_ts);PYLITH_CHECK_ERROR(err);

  TSEquationType eqType = TS_EQ_UNSPECIFIED;
  err = TSGetEquationType(_ts, &eqType);PYLITH_CHECK_ERROR(err);
  switch (eqType) {
  case TS_EQ_UNSPECIFIED: {
    throw std::logic_error("Unrecognized time stepping equation type for PETSc time stepping object.");
    break;
  } // unspecified
  case TS_EQ_EXPLCIIT:
  case TS_EQ_ODE_EXPLICIT:
  case TS_EQ_DAE_SEMI_EXPLICIT_INDEX1:
  case TS_EQ_DAE_SEMI_EXPLICIT_INDEX2:
  case TS_EQ_DAE_SEMI_EXPLICIT_INDEX3:
  case TS_EQ_DAE_SEMI_EXPLICIT_INDEHI:
    _formulationType = EXPLICIT;
    break;
  case TS_EQ_IMPLICIT:
  case TS_EQ_ODE_IMPLICIT:
  case TS_EQ_DAE_IMPLICIT_INDEX1:
  case TS_EQ_DAE_IMPLICIT_INDEX2:
  case TS_EQ_DAE_IMPLICIT_INDEX3:
  case TS_EQ_DAE_IMPLICIT_INDEXHI:
    _formultionType = IMPLICIT;
    break;
  default: {
  } // default
    assert(0);
    throw std::logic_error("Unknown PETSc time stepping equation type.");
  } // switch
    
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

  // Set initial solution.
  err = TSSetSolution(_ts, _solution->globalVector());PYLITH_CHECKE_ERROR(err);

  // Set callbacks.
  err = TSSetPreStep(_ts, prestep);PYLITH_CHECK_ERROR(err);
  err = TSSetPostStep(_ts, poststep);PYLITH_CHECK_ERROR(err);
  err = TSSetRHSJacobian(_ts, _jacobianRHS->matrix(), precondMatrix, reformRHSJacobian, (void*)this);PYLITH_CHECK_ERROR(err);
  err = TSSetRHSFunction(_ts, NULL, reformRHSResidual, (void*)this);PYLITH_CHECK_ERROR(err);
  
  if (IMPLICIT == _formulationType) {
    err = TSSetLHSFunction(_ts, NULL, reformLHSResidual, (void*)this);PYLITH_CHECK_ERROR(err);
    err = TSSetLHSJacobian(_ts, _jacobianLHS->matrix(), precondMatrix, reformLHSJacobian, (void*)this);PYLITH_CHECK_ERROR(err);
  } // if

  // Setup time stepper.
  err = TSSetUp(_ts);PYLITH_CHECK_ERROR(err);

  PYLITH_METHOD_END;
} // initialize

// ----------------------------------------------------------------------
// Solve time-dependent problem.
void
pylith::problems::TimeDependent::solve(void)
{ // solve
  PYLITH_METHOD_BEGIN;

  PetscErrorCode err = TSSolve(_ts, NULL);PYLITH_CHECK_ERROR(err);

  PYLITH_METHOD_END;
} // solve


// ----------------------------------------------------------------------
// Perform operations before advancing solution one time step.
void
pylith::problems::TimeDependent::prestep(void)
{ // prestep
  PYLITH_METHOD_BEGIN;

  // Get time and time step
  PetscErrorCode err;
  PylithReal dt;
  PylithReal t;
  err = TSGetTimeStep(_ts, &dt);PYLITH_CHECK_ERROR(err);
  err = TSGetTime(_ts, &t);

  // Set constraints.
  const size_t numConstraints = _constraints.size();
  for (size_t i=0; i < numConstraints; ++i) {
    _constraints[i]->setSolution(t, dt, _solution);
  } // for

  PYLITH_METHOD_END;
} // prestep

// ----------------------------------------------------------------------
// Perform operations after advancing solution one time step.
void
pylith::problems::TimeDependent::poststep(void)
{ // poststep
  PYLITH_METHOD_BEGIN;

  // :TODO: :INCOMPLETE:

  // Update state variables
  const size_t numIntegrators = _integrators.size();
  assert(numIntegrators > 0); // must have at least 1 integrator
  for (size_t i=0; i < numIntegrators; ++i) {
    _integrators[i]->updateStateVars(*_solution);
  } // for

  // Output

  PYLITH_METHOD_END;
} // poststep


// ----------------------------------------------------------------------
// Callback static method for reforming residual for RHS, G(t,u).
PetscErrorCode
pylith::problems::TimeDependent::reformRHSResidual(PetscTS ts,
						   PetscReal t,
						   PetscVec solutionVec,
						   PetscVec residualVec,
						   void* context)
{ // reformRHSResidual
  PYLITH_METHOD_BEGIN;

  // Get current time step.
  PylithReal dt;
  PetscErrorCode err = TSGetTimeStep(ts, &dt);PYLITH_CHECK_ERROR(err);

  pylith::problems::TimeDependent* problem = (pylith::problems::TimeDependent*)context;
  problem->reformRHSResidual(t, dt, solutionVec, residualVec);

  // If explicit time stepping, multiply RHS, G(t,u), by M^{-1}
  if (EXPLICIT == _formulationType) {
    assert(_jacobianLHS);

    // :KLUDGE: Should add check to see if we need to reform Jacobian
    problem->reformLHSJacobianExplicit(t, dt, solutionVec);

    PetscVec jacobianDiag = NULL;
    err = VecDuplicate(residualVec, &jacobianDiag);
    err = MatGetDiagonal(*_jacobianLHS, jacobianDiag);PYLITH_CHECK_ERROR(err);
    err = VecReciprocal(jacobianDiag);PYLITH_CHECK_ERROR(err);
    err = VecPointwiseMult(residualVec, jacobianDiag, residualVec);PYLITH_CHECK_ERROR(err);
  } // if

  PYLITH_METHOD_RETURN(0);
} // reformRHSResidual

  
// ----------------------------------------------------------------------
// Callback static method for reforming Jacobian for RHS, Jacobian of G(t,u).
PetscErrorCode
pylith::problems::TimeDependent::reformRHSJacobian(PetscTS ts,
						   PetscReal t,
						   PetscVec solutionVec,
						   PetscMat jacobianMat,
						   PetscMat precondMat,
						   void* context)
{ // reformRHSJacobian
  PYLITH_METHOD_BEGIN;

  // Get current time step.
  PylithReal dt;
  PetscErrorCode err = TSGetTimeStep(ts, &dt);PYLITH_CHECK_ERROR(err);

  pylith::problems::TimeDependent* problem = (pylith::problems::TimeDependent*)context;
  problem->reformRHSJacobian(t, dt, solutionVec, jacobianMat, precondMat);

  PYLITH_METHOD_RETURN(0);
} // reformRHSJacobian

// ----------------------------------------------------------------------
// Callback static method for reforming residual for LHS, F(t,u,\dot{u}).
PetscErrorCode
pylith::problems::TimeDependent::reformLHSResidual(PetscTS ts,
						   PetscReal t,
						   PetscVec solutionVec,
						   PetscVec residualVec,
						   void* context)
{ // reformLHSResidual
  PYLITH_METHOD_BEGIN;

  // Get current time step.
  PylithReal dt;
  PetscErrorCode err = TSGetTimeStep(ts, &dt);PYLITH_CHECK_ERROR(err);

  pylith::problems::TimeDependent* problem = (pylith::problems::TimeDependent*)context;
  problem->reformLHSResidual(t, dt, solutionVec, residualVec);

  PYLITH_METHOD_RETURN(0);
} // reformLHSResidual

  
// ----------------------------------------------------------------------
// Callback static method for reforming Jacobian for LHS, Jacobian of F(t,u,\dot{u}).
PetscErrorCode
pylith::problems::TimeDependent::reformLHSJacobian(PetscTS ts,
						   PetscReal t,
						   PetscVec solutionVec,
						   PetscMat jacobianMat,
						   PetscMat precondMat,
						   void* context)
{ // reformLHSJacobian
  PYLITH_METHOD_BEGIN;

  // Get current time step.
  PylithReal dt;
  PetscErrorCode err = TSGetTimeStep(ts, &dt);PYLITH_CHECK_ERROR(err);

  pylith::problems::TimeDependent* problem = (pylith::problems::TimeDependent*)context;
  problem->reformLHSJacobianImplicit(t, dt, solutionVec, jacobianMat, precondMat);

  PYLITH_METHOD_RETURN(0);
} // reformLHSJacobian


// ----------------------------------------------------------------------
// Callback static method for operations before advancing solution one time step.
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
// Callback static method for operations after advancing solution one time step.
void
pylith::problems::TimeDependent::poststep(PetscTS ts)
{ // poststep
  PYLITH_METHOD_BEGIN;

  TimeDependent* problem = NULL;
  PetscErrorCode err = TSGetApplicationContext(ts, (void*)problem);PYLITH_CHECK_ERROR(err);assert(problem);
  problem->poststep();

  PYLITH_METHOD_RETURN(0);
} // poststep


// End of file
