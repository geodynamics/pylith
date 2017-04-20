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
#include "pylith/feassemble/ConstraintPointwise.hh" // USES ConstraintPointwise

#include "petscts.h" // USES PetscTS

#include "pylith/utils/error.hh" // USES PYLITH_CHECK_ERROR
#include "pylith/utils/journals.hh" // USES PYLITH_JOURNAL_*
#include <cassert> // USES assert()

// ----------------------------------------------------------------------
const char* pylith::problems::TimeDependent::_pyreComponent = "timedependent";

// ----------------------------------------------------------------------
// Constructor
pylith::problems::TimeDependent::TimeDependent(void) :
    _startTime(0.0),
    _dtInitial(1.0),
    _totalTime(0.0),
    _maxTimeSteps(0),
    _ts(0),
    _formulationType(IMPLICIT)
{ // constructor
    PyreComponent::name(_pyreComponent);
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

    Problem::deallocate();
    
    PetscErrorCode err = TSDestroy(&_ts); PYLITH_CHECK_ERROR(err);

    PYLITH_METHOD_END;
} // deallocate

// ----------------------------------------------------------------------
// Set start time for problem.
void
pylith::problems::TimeDependent::startTime(const double value)
{ // startTime
    PYLITH_JOURNAL_DEBUG("startTime(value="<<value<<")");

    _startTime = value;
} // startTime

// ----------------------------------------------------------------------
// Get start time for problem.
double
pylith::problems::TimeDependent::startTime(void) const
{ // startTime
    return _startTime;
} // startTime

// ----------------------------------------------------------------------
// Set total time for problem.
void
pylith::problems::TimeDependent::totalTime(const double value)
{ // totalTime
    PYLITH_METHOD_BEGIN;
    PYLITH_JOURNAL_DEBUG("totalTime(value="<<value<<")");

    if (value < 0.0) {
        std::ostringstream msg;
        msg << "Total time (nondimensional) for problem (" << value << ") must be positive.";
        throw std::runtime_error(msg.str());
    }     // if
    _totalTime = value;

    PYLITH_METHOD_END;
} // totalTime

// ----------------------------------------------------------------------
// Get total time for problem.
double
pylith::problems::TimeDependent::totalTime(void) const
{ // totalTime
    return _totalTime;
} // totalTime

// ----------------------------------------------------------------------
// Set maximum number of time steps.
void
pylith::problems::TimeDependent::maxTimeSteps(const size_t value)
{ // maxTimeSteps
    PYLITH_METHOD_BEGIN;
    PYLITH_JOURNAL_DEBUG("maxTimeSteps(value="<<value<<")");

    if (value <= 0) {
        std::ostringstream msg;
        msg << "Maximum number of time teps for problem (" << value << ") must be positive.";
        throw std::runtime_error(msg.str());
    }     // if
    _maxTimeSteps = value;

    PYLITH_METHOD_END;
} // maxTimeSteps

// ----------------------------------------------------------------------
// Get maximum number of time steps.
size_t
pylith::problems::TimeDependent::maxTimeSteps(void) const
{ // maxTimeSteps
    return _maxTimeSteps;
} // maxTimeSteps

// ----------------------------------------------------------------------
// Set initial time step for problem.
void
pylith::problems::TimeDependent::dtInitial(const double value)
{ // dtInitial
    PYLITH_METHOD_BEGIN;
    PYLITH_JOURNAL_DEBUG("dtInitial(value="<<value<<")");

    if (value < 0.0) {
        std::ostringstream msg;
        msg << "Initial time step (nondimensional) for problem (" << value << ") must be positive.";
        throw std::runtime_error(msg.str());
    }     // if
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
pylith::problems::TimeDependent::initialize(void)
{ // initialize
    PYLITH_METHOD_BEGIN;
    PYLITH_JOURNAL_DEBUG("initialize()");

    Problem::initialize();

    assert(_solution);

    PetscErrorCode err = TSDestroy(&_ts); PYLITH_CHECK_ERROR(err); assert(!_ts);
    const pylith::topology::Mesh& mesh = _solution->mesh();
    err = TSCreate(mesh.comm(), &_ts); PYLITH_CHECK_ERROR(err); assert(_ts);
    err = TSSetType(_ts, TSBEULER); PYLITH_CHECK_ERROR(err); // Backward Euler is default time stepping method.
    err = TSSetExactFinalTime(_ts, TS_EXACTFINALTIME_STEPOVER); PYLITH_CHECK_ERROR(err); // Ok to step over final time.
    err = TSSetFromOptions(_ts); PYLITH_CHECK_ERROR(err);
    err = TSSetApplicationContext(_ts, (void*)this); PYLITH_CHECK_ERROR(err);

    #if 0
    TSEquationType eqType = TS_EQ_UNSPECIFIED;
    err = TSGetEquationType(_ts, &eqType); PYLITH_CHECK_ERROR(err);
    switch (eqType) {
    case TS_EQ_UNSPECIFIED: {
        PYLITH_JOURNAL_ERROR("Unspecified time stepping equation type for PETSc time stepping object.");
        throw std::logic_error("Unspecified time stepping equation type for PETSc time stepping object.");
        break;
    }         // unspecified
    case TS_EQ_EXPLICIT:
    case TS_EQ_ODE_EXPLICIT:
    case TS_EQ_DAE_SEMI_EXPLICIT_INDEX1:
    case TS_EQ_DAE_SEMI_EXPLICIT_INDEX2:
    case TS_EQ_DAE_SEMI_EXPLICIT_INDEX3:
    case TS_EQ_DAE_SEMI_EXPLICIT_INDEXHI:
        PYLITH_JOURNAL_DEBUG("Recognized PetscTS as explicit time stepping.");
        _formulationType = EXPLICIT;
        break;
    case TS_EQ_IMPLICIT:
    case TS_EQ_ODE_IMPLICIT:
    case TS_EQ_DAE_IMPLICIT_INDEX1:
    case TS_EQ_DAE_IMPLICIT_INDEX2:
    case TS_EQ_DAE_IMPLICIT_INDEX3:
    case TS_EQ_DAE_IMPLICIT_INDEXHI:
        PYLITH_JOURNAL_DEBUG("Recognized PetscTS as implicit time stepping.");
        _formulationType = IMPLICIT;
        break;
    default:
        PYLITH_JOURNAL_DEBUG("Unknown PETSc time stepping equation type.");
        throw std::logic_error("Unknown PETSc time stepping equation type.");
    }         // switch
        #else
    TSType tsType = TSBEULER;
    err = TSGetType(_ts, &tsType); PYLITH_CHECK_ERROR(err);
    std::string tsString = tsType;
    if (tsString == std::string(TSEULER)
        || tsString == std::string(TSSSP)
        || tsString == std::string(TSRK)) {
        PYLITH_JOURNAL_DEBUG("Recognized PetscTS as explicit time stepping.");
        _formulationType = EXPLICIT;
    } else if (tsString == std::string(TSBEULER)
               || tsString == std::string(TSCN)
               || tsString == std::string(TSTHETA)
               || tsString == std::string(TSALPHA)
               || tsString == std::string(TSROSW)) {
        PYLITH_JOURNAL_DEBUG("Recognized PetscTS as implicit time stepping.");
        _formulationType = IMPLICIT;
    } else {
        // TSPSEUDO:
        // TSSUNDIALS:
        // TSPYTHON:
        // TSALPHA2:
        // TSGLLE:
        // TSGLEE:
        // TSARKIMEX:
        // TSEIMEX:
        // TSMIMEX:
        // TSBDF:
        PYLITH_JOURNAL_ERROR("Unable to determine if PETSc TS type '"<<tsType<<"' is implicit or explicit.");
        throw std::logic_error("Unable to determine if PETSc TS type is implicit or explicit.");
    }         // if/else
        #endif

    // Set time stepping paramters.
    switch (this->solverType()) {
    case LINEAR:
        PYLITH_JOURNAL_DEBUG("Setting PetscTS problem type to 'linear'.");
        err = TSSetProblemType(_ts, TS_LINEAR); PYLITH_CHECK_ERROR(err);
        break;
    case NONLINEAR:
        PYLITH_JOURNAL_DEBUG("Setting PetscTS problem type to 'nonlinear'.");
        err = TSSetProblemType(_ts, TS_NONLINEAR); PYLITH_CHECK_ERROR(err);
        break;
    default:
        PYLITH_JOURNAL_ERROR("Unknown problem type.");
        throw std::logic_error("Unknown problem type.");
    }     // switch
    PYLITH_JOURNAL_DEBUG("Setting PetscTS parameters: dtInitial="<<_dtInitial<<", startTime="<<_startTime<<", maxTimeSteps="<<_maxTimeSteps<<", totalTime="<<_totalTime);
    err = TSSetInitialTimeStep(_ts, _startTime, _dtInitial); PYLITH_CHECK_ERROR(err);
    err = TSSetDuration(_ts, _maxTimeSteps, _totalTime); PYLITH_CHECK_ERROR(err);

    // Set initial solution.
    _solution->zeroLocal();
    PYLITH_JOURNAL_ERROR(":TODO: @brad Implement setting initial solution.");
    // :TODO: Set initial conditions.
    PetscVec solutionVec;
    err = DMCreateGlobalVector(_solution->dmMesh(), &solutionVec); PYLITH_CHECK_ERROR(err);
    _solution->scatterLocalToVector(solutionVec);
    PYLITH_JOURNAL_DEBUG("Setting PetscTS initial conditions using global vector for solution.");
    err = TSSetSolution(_ts, solutionVec); PYLITH_CHECK_ERROR(err);
    err = VecDestroy(&solutionVec); PYLITH_CHECK_ERROR(err);

    // Set callbacks.
    PYLITH_JOURNAL_DEBUG("Setting PetscTS callbacks prestep(), poststep(), computeRHSJacobian(), and computeRHSFunction().");
    err = TSSetPreStep(_ts, prestep); PYLITH_CHECK_ERROR(err);
    err = TSSetPostStep(_ts, poststep); PYLITH_CHECK_ERROR(err);
    err = TSSetRHSJacobian(_ts, NULL, NULL, computeRHSJacobian, (void*)this); PYLITH_CHECK_ERROR(err);
    err = TSSetRHSFunction(_ts, NULL, computeRHSResidual, (void*)this); PYLITH_CHECK_ERROR(err);

    if (IMPLICIT == _formulationType) {
        PYLITH_JOURNAL_DEBUG("Setting PetscTS callbacks computeLHSJacobian(), and computeLHSFunction().");
        err = TSSetIFunction(_ts, NULL, computeLHSResidual, (void*)this); PYLITH_CHECK_ERROR(err);
        err = TSSetIJacobian(_ts, NULL, NULL, computeLHSJacobian, (void*)this); PYLITH_CHECK_ERROR(err);
    } // if

    // Setup time stepper.
    err = TSSetUp(_ts); PYLITH_CHECK_ERROR(err);

    // Setup field to hold inverse of lumped LHS Jacobian (if explicit).
    if (EXPLICIT == _formulationType) {
        PYLITH_JOURNAL_DEBUG("Setting up field for inverse of lumped LHS Jacobian.");

        delete _jacobianLHSLumpedInv; _jacobianLHSLumpedInv = new pylith::topology::Field(_solution->mesh()); assert(_jacobianLHSLumpedInv);
        _jacobianLHSLumpedInv->cloneSection(*_solution);
    } // if

    PYLITH_METHOD_END;
} // initialize

// ----------------------------------------------------------------------
// Solve time-dependent problem.
void
pylith::problems::TimeDependent::solve(void)
{ // solve
    PYLITH_METHOD_BEGIN;
    PYLITH_JOURNAL_DEBUG("solve()");

    PetscErrorCode err = TSSolve(_ts, NULL); PYLITH_CHECK_ERROR(err);

    PYLITH_METHOD_END;
} // solve


// ----------------------------------------------------------------------
// Perform operations before advancing solution one time step.
void
pylith::problems::TimeDependent::prestep(void)
{ // prestep
    PYLITH_METHOD_BEGIN;
    PYLITH_JOURNAL_DEBUG("prestep()");

    // Get time and time step
    PetscErrorCode err;
    PylithReal dt;
    PylithReal t;
    err = TSGetTimeStep(_ts, &dt); PYLITH_CHECK_ERROR(err);
    err = TSGetTime(_ts, &t); PYLITH_CHECK_ERROR(err);

    // Prepare constraints for a new time step.
    const size_t numConstraints = _constraints.size();
    for (size_t i=0; i < numConstraints; ++i) {
        _constraints[i]->prestep(t, dt);
    } // for

    // Prepare integrators for a new time step.
    const size_t numIntegrators = _integrators.size();
    for (size_t i=0; i < numIntegrators; ++i) {
        _integrators[i]->prestep(t, dt);
    } // for

    PYLITH_METHOD_END;
} // prestep

// ----------------------------------------------------------------------
// Perform operations after advancing solution one time step.
void
pylith::problems::TimeDependent::poststep(void)
{ // poststep
    PYLITH_METHOD_BEGIN;
    PYLITH_JOURNAL_DEBUG("poststep()");

    // Get current solution.
    // :QUESTION: :MATT: What time does this solution correspond to?
    PetscErrorCode err;
    PylithReal t;
    PylithInt dt;
    PetscVec solutionVec = NULL;
    err = TSGetTime(_ts, &t); PYLITH_CHECK_ERROR(err);
    err = TSGetTotalSteps(_ts, &dt); PYLITH_CHECK_ERROR(err);
    err = TSGetSolution(_ts, &solutionVec); PYLITH_CHECK_ERROR(err);

    // Update PyLith view of the solution.
    assert(_solution);
    _solution->scatterVectorToLocal(solutionVec);

    // Output solution.
    const size_t numOutput = _outputs.size();
    for (size_t i=0; i < numOutput; ++i) {
#if 1
        PYLITH_JOURNAL_ERROR(":TODO: @brad Implement solution output in poststep().");
#else
        _outputs[i].writeTimeStep(t, dt, *_solution);
#endif
    } // for

    // Update state variables
    const size_t numIntegrators = _integrators.size();
    assert(numIntegrators > 0);     // must have at least 1 integrator
    for (size_t i=0; i < numIntegrators; ++i) {
        _integrators[i]->updateStateVars(*_solution);
#if 1
        PYLITH_JOURNAL_ERROR(":TODO: @brad Implement integrator output in poststep().");
#else
        _integrators[i]->writeTimeStep(t, *_solution);
#endif
    }  // for

    PYLITH_METHOD_END;
} // poststep


// ----------------------------------------------------------------------
// Callback static method for computing residual for RHS, G(t,s).
PetscErrorCode
pylith::problems::TimeDependent::computeRHSResidual(PetscTS ts,
                                                    PetscReal t,
                                                    PetscVec solutionVec,
                                                    PetscVec residualVec,
                                                    void* context)
{ // computeRHSResidual
    PYLITH_METHOD_BEGIN;

    journal::debug_t debug(_pyreComponent);
    debug << journal::at(__HERE__)
          << "computeRHSResidual(ts="<<ts<<", t="<<t<<", solutionVec="<<solutionVec<<", residualVec="<<residualVec<<", context="<<context<<")" << journal::endl;

    // Get current time step.
    PylithReal dt;
    PetscErrorCode err = TSGetTimeStep(ts, &dt); PYLITH_CHECK_ERROR(err);

    pylith::problems::TimeDependent* problem = (pylith::problems::TimeDependent*)context;
    problem->Problem::computeRHSResidual(residualVec, t, dt, solutionVec);

    // If explicit time stepping, multiply RHS, G(t,s), by M^{-1}
    if (EXPLICIT == problem->_formulationType) {
        problem->Problem::computeLHSJacobianLumpedInv(t, dt, solutionVec);

        assert(problem->_jacobianLHSLumpedInv);
        err = VecPointwiseMult(residualVec, problem->_jacobianLHSLumpedInv->localVector(), residualVec); PYLITH_CHECK_ERROR(err);
    }     // if

    PYLITH_METHOD_RETURN(0);
} // computeRHSResidual


// ----------------------------------------------------------------------
// Callback static method for computing Jacobian for RHS, Jacobian of G(t,s).
PetscErrorCode
pylith::problems::TimeDependent::computeRHSJacobian(PetscTS ts,
                                                    PetscReal t,
                                                    PetscVec solutionVec,
                                                    PetscMat jacobianMat,
                                                    PetscMat precondMat,
                                                    void* context)
{ // computeRHSJacobian
    PYLITH_METHOD_BEGIN;

    journal::debug_t debug(_pyreComponent);
    debug << journal::at(__HERE__)
          << "computeRHSJacobian(ts="<<ts<<", t="<<t<<", solutionVec="<<solutionVec<<", jacobianMat="<<jacobianMat<<", precondMat="<<precondMat<<", context="<<context<<")" << journal::endl;

    // Get current time step.
    PylithReal dt;
    PetscErrorCode err = TSGetTimeStep(ts, &dt); PYLITH_CHECK_ERROR(err);

    pylith::problems::TimeDependent* problem = (pylith::problems::TimeDependent*)context;
    problem->Problem::computeRHSJacobian(jacobianMat, precondMat, t, dt, solutionVec);

    PYLITH_METHOD_RETURN(0);
} // computeRHSJacobian

// ----------------------------------------------------------------------
// Callback static method for computing residual for LHS, F(t,s,\dot{s}).
PetscErrorCode
pylith::problems::TimeDependent::computeLHSResidual(PetscTS ts,
                                                    PetscReal t,
                                                    PetscVec solutionVec,
                                                    PetscVec solutionDotVec,
                                                    PetscVec residualVec,
                                                    void* context)
{ // computeLHSResidual
    PYLITH_METHOD_BEGIN;

    journal::debug_t debug(_pyreComponent);
    debug << journal::at(__HERE__)
          << "computeLHSResidual(ts="<<ts<<", t="<<t<<", solutionVec="<<solutionVec<<", solutionDotVec="<<solutionDotVec<<", residualVec="<<residualVec<<", context="<<context<<")" << journal::endl;

    // Get current time step.
    PylithReal dt;
    PetscErrorCode err = TSGetTimeStep(ts, &dt); PYLITH_CHECK_ERROR(err);

    pylith::problems::TimeDependent* problem = (pylith::problems::TimeDependent*)context;
    problem->Problem::computeLHSResidual(residualVec, t, dt, solutionVec, solutionDotVec);

    PYLITH_METHOD_RETURN(0);
} // computeLHSResidual


// ----------------------------------------------------------------------
// Callback static method for computing Jacobian for LHS, Jacobian of F(t,s,\dot{s}).
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

    journal::debug_t debug(_pyreComponent);
    debug << journal::at(__HERE__)
          << "computeLHSJacobian(ts="<<ts<<", t="<<t<<", solutionVec="<<solutionVec<<", solutionDotVec="<<solutionDotVec<<", tshift="<<tshift<<", jacobianMat="<<jacobianMat<<", precondMat="<<precondMat<<", context="<<context<<")" <<
    journal::endl;

    // Get current time step.
    PylithReal dt;
    PetscErrorCode err = TSGetTimeStep(ts, &dt); PYLITH_CHECK_ERROR(err);

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

    journal::debug_t debug(_pyreComponent);
    debug << journal::at(__HERE__)
          << "prestep(ts="<<ts<<")" << journal::endl;

    TimeDependent* problem = NULL;
    PetscErrorCode err = TSGetApplicationContext(ts, (void*)&problem); PYLITH_CHECK_ERROR(err); assert(problem);
    problem->prestep();

    PYLITH_METHOD_RETURN(0);
} // prestep


// ----------------------------------------------------------------------
// Callback static method for operations after advancing solution one time step.
PetscErrorCode
pylith::problems::TimeDependent::poststep(PetscTS ts)
{ // poststep
    PYLITH_METHOD_BEGIN;

    journal::debug_t debug(_pyreComponent);
    debug << journal::at(__HERE__)
          << "poststep(ts="<<ts<<")" << journal::endl;

    TimeDependent* problem = NULL;
    PetscErrorCode err = TSGetApplicationContext(ts, (void*)&problem); PYLITH_CHECK_ERROR(err); assert(problem);
    problem->poststep();

    PYLITH_METHOD_RETURN(0);
} // poststep


// End of file
