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

#include "pylith/topology/Mesh.hh" // USES Mesh
#include "pylith/topology/Field.hh" // USES Field

#include "pylith/feassemble/IntegratorPointwise.hh" // USES IntegratorPointwise
#include "pylith/feassemble/Constraint.hh" // USES Constraint
#include "pylith/topology/Jacobian.hh" // USES Jacobian

#include "petscts.h" // USES PetscTS

#include "pylith/utils/error.hh" // USES PYLITH_CHECK_ERROR
#include "journal/debug.h" // USES journal::debug_t
#include "journal/error.h" // USES journal::error_t
#include <cassert> // USES assert()

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

  PetscErrorCode err = TSDestroy(&_ts);PYLITH_CHECK_ERROR(err);

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
// Set maximum number of time steps.
void
pylith::problems::TimeDependent::maxTimeSteps(const PetscInt value)
{ // maxTimeSteps
  PYLITH_METHOD_BEGIN;

  if (value <= 0) {
    std::ostringstream msg;
    msg << "Maximum number of time teps for problem (" << value << ") must be positive.";
    throw std::runtime_error(msg.str());
  } // if
  _maxTimeSteps = value;

  PYLITH_METHOD_END;  
} // maxTimeSteps

// ----------------------------------------------------------------------
// Get maximum number of time steps.
PetscInt
pylith::problems::TimeDependent::maxTimeSteps(void) const
{ // maxTimeSteps
  return _maxTimeSteps;
} // maxTimeSteps

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
pylith::problems::TimeDependent::initialize(pylith::topology::Field* solution)
{ // initialize
  PYLITH_METHOD_BEGIN;

  journal::debug_t debug("problem");
  debug << journal::at(__HERE__)
	<< "TimeDependent::initialize(solution="<<solution<< ")" << journal::endl;

  assert(solution);

  _solution = solution;

  PetscErrorCode err = TSDestroy(&_ts);PYLITH_CHECK_ERROR(err);assert(!_ts);
  const pylith::topology::Mesh& mesh = solution->mesh();
  err = TSCreate(mesh.comm(), &_ts);PYLITH_CHECK_ERROR(err);assert(_ts);
  err = TSSetFromOptions(_ts);PYLITH_CHECK_ERROR(err);

  TSEquationType eqType = TS_EQ_UNSPECIFIED;
  err = TSGetEquationType(_ts, &eqType);PYLITH_CHECK_ERROR(err);
  switch (eqType) {
  case TS_EQ_UNSPECIFIED: {
    throw std::logic_error("Unrecognized time stepping equation type for PETSc time stepping object.");
    break;
  } // unspecified
  case TS_EQ_EXPLICIT:
  case TS_EQ_ODE_EXPLICIT:
  case TS_EQ_DAE_SEMI_EXPLICIT_INDEX1:
  case TS_EQ_DAE_SEMI_EXPLICIT_INDEX2:
  case TS_EQ_DAE_SEMI_EXPLICIT_INDEX3:
  case TS_EQ_DAE_SEMI_EXPLICIT_INDEXHI:
    _formulationType = EXPLICIT;
    break;
  case TS_EQ_IMPLICIT:
  case TS_EQ_ODE_IMPLICIT:
  case TS_EQ_DAE_IMPLICIT_INDEX1:
  case TS_EQ_DAE_IMPLICIT_INDEX2:
  case TS_EQ_DAE_IMPLICIT_INDEX3:
  case TS_EQ_DAE_IMPLICIT_INDEXHI:
    _formulationType = IMPLICIT;
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
  err = TSSetInitialTimeStep(_ts, _dtInitial, _startTime);PYLITH_CHECK_ERROR(err);
  err = TSSetDuration(_ts, _maxTimeSteps, _totalTime);PYLITH_CHECK_ERROR(err);

  // Set initial solution.
  err = TSSetSolution(_ts, _solution->globalVector());PYLITH_CHECK_ERROR(err);

  // Set callbacks.
  err = TSSetPreStep(_ts, prestep);PYLITH_CHECK_ERROR(err);
  err = TSSetPostStep(_ts, poststep);PYLITH_CHECK_ERROR(err);
  err = TSSetRHSJacobian(_ts, NULL, NULL, computeRHSJacobian, (void*)this);PYLITH_CHECK_ERROR(err);
  err = TSSetRHSFunction(_ts, NULL, computeRHSResidual, (void*)this);PYLITH_CHECK_ERROR(err);
  
  if (IMPLICIT == _formulationType) {
    err = TSSetIFunction(_ts, NULL, computeLHSResidual, (void*)this);PYLITH_CHECK_ERROR(err);
    err = TSSetIJacobian(_ts, NULL, NULL, computeLHSJacobian, (void*)this);PYLITH_CHECK_ERROR(err);
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

  journal::debug_t debug("problem");
  debug << journal::at(__HERE__)
	<< "TimeDependent::solve()" << journal::endl;

  PetscErrorCode err = TSSolve(_ts, NULL);PYLITH_CHECK_ERROR(err);

  PYLITH_METHOD_END;
} // solve


// ----------------------------------------------------------------------
// Perform operations before advancing solution one time step.
void
pylith::problems::TimeDependent::prestep(void)
{ // prestep
  PYLITH_METHOD_BEGIN;

  journal::debug_t debug("problem");
  debug << journal::at(__HERE__)
	<< "TimeDependent::prestep()" << journal::endl;

  // Get time and time step
  PetscErrorCode err;
  PylithReal dt;
  PylithReal t;
  err = TSGetTimeStep(_ts, &dt);PYLITH_CHECK_ERROR(err);
  err = TSGetTime(_ts, &t);

  journal::error_t error("problem");
  error << journal::at(__HERE__)
	  << "TimeDependent::prestep() missing setting of constraints; "
	  << "Constraint object missing setAuxFields(t, dt)." << journal::endl;

  // Set constraints.
#if 0 // :KLUDGE: :TODO: Implement this.
  const size_t numConstraints = _constraints.size();
  for (size_t i=0; i < numConstraints; ++i) {
    _constraints[i]->setAuxFields(t, dt);
  } // for
#endif

  PYLITH_METHOD_END;
} // prestep

// ----------------------------------------------------------------------
// Perform operations after advancing solution one time step.
void
pylith::problems::TimeDependent::poststep(void)
{ // poststep
  PYLITH_METHOD_BEGIN;

  journal::debug_t debug("problem");
  debug << journal::at(__HERE__)
	<< "TimeDependent::poststep()" << journal::endl;
  
  // Get current solution.
  // :QUESTION: :MATT: What time does this solution correspond to?
  PetscVec solutionVec = NULL;
  PetscErrorCode err = TSGetSolution(_ts, &solutionVec);PYLITH_CHECK_ERROR(err);

  // Update PyLith view of the solution.
  assert(_solution);
  _solution->scatterGlobalToLocal(solutionVec);

  // Update state variables
  const size_t numIntegrators = _integrators.size();
  assert(numIntegrators > 0); // must have at least 1 integrator
  for (size_t i=0; i < numIntegrators; ++i) {
    _integrators[i]->updateStateVars(*_solution);
  } // for

  journal::error_t error("problem");
  error << journal::at(__HERE__)
	  << "TimeDependent::poststep() missing output." << journal::endl;

  // :TODO: Output [this is output at whatever time the solution corresponds to.]

  PYLITH_METHOD_END;
} // poststep


// ----------------------------------------------------------------------
// Callback static method for computeing residual for RHS, G(t,s).
PetscErrorCode
pylith::problems::TimeDependent::computeRHSResidual(PetscTS ts,
						    PetscReal t,
						    PetscVec solutionVec,
						    PetscVec residualVec,
						    void* context)
{ // computeRHSResidual
  PYLITH_METHOD_BEGIN;

  journal::debug_t debug("problem");
  debug << journal::at(__HERE__)
	<< "TimeDependent::computeRHSResidual(ts="<<ts<<", t="<<t<<", solutionVec="<<solutionVec<<", residualVec="<<residualVec<<", context="<<context<<")" << journal::endl;
  
  // Get current time step.
  PylithReal dt;
  PetscErrorCode err = TSGetTimeStep(ts, &dt);PYLITH_CHECK_ERROR(err);

  pylith::problems::TimeDependent* problem = (pylith::problems::TimeDependent*)context;
  problem->Problem::computeRHSResidual(residualVec, t, dt, solutionVec);

  // If explicit time stepping, multiply RHS, G(t,s), by M^{-1}
  if (EXPLICIT == problem->_formulationType) {

    // :KLUDGE: :TODO: Should add check to see if we need to compute Jacobian
    problem->Problem::computeLHSJacobianLumpedInv(t, dt, solutionVec);

    assert(problem->_jacobianLHSLumpedInv);
    err = VecPointwiseMult(residualVec, problem->_jacobianLHSLumpedInv->localVector(), residualVec);PYLITH_CHECK_ERROR(err);
  } // if

  PYLITH_METHOD_RETURN(0);
} // computeRHSResidual

  
// ----------------------------------------------------------------------
// Callback static method for computeing Jacobian for RHS, Jacobian of G(t,s).
PetscErrorCode
pylith::problems::TimeDependent::computeRHSJacobian(PetscTS ts,
						    PetscReal t,
						    PetscVec solutionVec,
						    PetscMat jacobianMat,
						    PetscMat precondMat,
						    void* context)
{ // computeRHSJacobian
  PYLITH_METHOD_BEGIN;

  journal::debug_t debug("problem");
  debug << journal::at(__HERE__)
	<< "TimeDependent::computeRHSJacobian(ts="<<ts<<", t="<<t<<", solutionVec="<<solutionVec<<", jacobianMat="<<jacobianMat<<", precondMat="<<precondMat<<", context="<<context<<")" << journal::endl;

  // Get current time step.
  PylithReal dt;
  PetscErrorCode err = TSGetTimeStep(ts, &dt);PYLITH_CHECK_ERROR(err);

  pylith::problems::TimeDependent* problem = (pylith::problems::TimeDependent*)context;
  problem->Problem::computeRHSJacobian(jacobianMat, precondMat, t, dt, solutionVec);

  PYLITH_METHOD_RETURN(0);
} // computeRHSJacobian

// ----------------------------------------------------------------------
// Callback static method for computeing residual for LHS, F(t,s,\dot{s}).
PetscErrorCode
pylith::problems::TimeDependent::computeLHSResidual(PetscTS ts,
						    PetscReal t,
						    PetscVec solutionVec,
						    PetscVec solutionDotVec,
						    PetscVec residualVec,
						    void* context)
{ // computeLHSResidual
  PYLITH_METHOD_BEGIN;

  journal::debug_t debug("problem");
  debug << journal::at(__HERE__)
	<< "TimeDependent::computeLHSResidual(ts="<<ts<<", t="<<t<<", solutionVec="<<solutionVec<<", solutionDotVec="<<solutionDotVec<<", residualVec="<<residualVec<<", context="<<context<<")" << journal::endl;
  
  // Get current time step.
  PylithReal dt;
  PetscErrorCode err = TSGetTimeStep(ts, &dt);PYLITH_CHECK_ERROR(err);

  pylith::problems::TimeDependent* problem = (pylith::problems::TimeDependent*)context;
  problem->Problem::computeLHSResidual(residualVec, t, dt, solutionVec, solutionDotVec);

  PYLITH_METHOD_RETURN(0);
} // computeLHSResidual

  
// ----------------------------------------------------------------------
// Callback static method for computeing Jacobian for LHS, Jacobian of F(t,s,\dot{s}).
PetscErrorCode
pylith::problems::TimeDependent::computeLHSJacobian(PetscTS ts,
						    PetscReal t,
						    PetscVec solutionVec,
						    PetscVec solutionDotVec,
						    PetscReal tshift,
						    PetscMat jacobianMat,
						    PetscMat precondMat,
						    void* context)
{ // computeLHSJacobian
  PYLITH_METHOD_BEGIN;

  journal::debug_t debug("problem");
  debug << journal::at(__HERE__)
	<< "TimeDependent::computeLHSJacobian(ts="<<ts<<", t="<<t<<", solutionVec="<<solutionVec<<", solutionDotVec="<<solutionDotVec<<", tshift="<<tshift<<", jacobianMat="<<jacobianMat<<", precondMat="<<precondMat<<", context="<<context<<")" << journal::endl;

  // Get current time step.
  PylithReal dt;
  PetscErrorCode err = TSGetTimeStep(ts, &dt);PYLITH_CHECK_ERROR(err);

  pylith::problems::TimeDependent* problem = (pylith::problems::TimeDependent*)context;
  problem->computeLHSJacobianImplicit(jacobianMat, precondMat, t, dt, tshift, solutionVec, solutionDotVec);

  PYLITH_METHOD_RETURN(0);
} // computeLHSJacobian


// ----------------------------------------------------------------------
// Callback static method for operations before advancing solution one time step.
PetscErrorCode
pylith::problems::TimeDependent::prestep(PetscTS ts)
{ // prestep
  PYLITH_METHOD_BEGIN;

  journal::debug_t debug("problem");
  debug << journal::at(__HERE__)
	<< "TimeDependent::prestep(ts="<<ts<<")" << journal::endl;

  TimeDependent* problem = NULL;
  PetscErrorCode err = TSGetApplicationContext(ts, (void*)problem);PYLITH_CHECK_ERROR(err);assert(problem);
  problem->prestep();

  PYLITH_METHOD_RETURN(0);
} // prestep


// ----------------------------------------------------------------------
// Callback static method for operations after advancing solution one time step.
PetscErrorCode
pylith::problems::TimeDependent::poststep(PetscTS ts)
{ // poststep
  PYLITH_METHOD_BEGIN;

  journal::debug_t debug("problem");
  debug << journal::at(__HERE__)
	<< "TimeDependent::poststep(ts="<<ts<<")" << journal::endl;

  TimeDependent* problem = NULL;
  PetscErrorCode err = TSGetApplicationContext(ts, (void*)problem);PYLITH_CHECK_ERROR(err);assert(problem);
  problem->poststep();

  PYLITH_METHOD_RETURN(0);
} // poststep


// End of file
