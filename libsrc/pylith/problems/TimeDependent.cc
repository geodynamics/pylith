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

#include "pylith/feassemble/Integrator.hh" // USES Integrator
#include "pylith/feassemble/Constraint.hh" // USES Constraint
#include "pylith/problems/ObserversSoln.hh" // USES ObserversSoln
#include "pylith/problems/InitialCondition.hh" // USES InitialCondition
#include "pylith/problems/ProgressMonitorTime.hh" // USES ProgressMonitorTime

#include "spatialdata/units/Nondimensional.hh" // USES Nondimensional

#include "petscts.h" // USES PetscTS

#include "pylith/utils/error.hh" // USES PYLITH_CHECK_ERROR
#include "pylith/utils/journals.hh" // USES PYLITH_COMPONENT_*
#include <cassert> // USES assert()

// ---------------------------------------------------------------------------------------------------------------------
namespace pylith {
    namespace problems {
        class _TimeDependent {
public:

            static const char* pyreComponent;
        }; // _TimeDependent

        const char* _TimeDependent::pyreComponent = "timedependent";
    } // problems
} // pylith

// ---------------------------------------------------------------------------------------------------------------------
// Constructor
pylith::problems::TimeDependent::TimeDependent(void) :
    _startTime(0.0),
    _endTime(0.0),
    _dtInitial(1.0),
    _maxTimeSteps(0),
    _ts(NULL),
    _monitor(NULL),
    _shouldNotifyIC(false) {
    PyreComponent::setName(_TimeDependent::pyreComponent);
} // constructor


// ---------------------------------------------------------------------------------------------------------------------
// Destructor
pylith::problems::TimeDependent::~TimeDependent(void) {
    deallocate();
} // destructor


// ---------------------------------------------------------------------------------------------------------------------
// Deallocate PETSc and local data structures.
void
pylith::problems::TimeDependent::deallocate(void) {
    PYLITH_METHOD_BEGIN;

    Problem::deallocate();

    PetscErrorCode err = TSDestroy(&_ts);PYLITH_CHECK_ERROR(err);

    _monitor = NULL; // Memory handle in Python. :TODO: Use shared pointer.

    PYLITH_METHOD_END;
} // deallocate


// ---------------------------------------------------------------------------------------------------------------------
// Set start time for problem.
void
pylith::problems::TimeDependent::setStartTime(const double value) {
    PYLITH_COMPONENT_DEBUG("startTime(value="<<value<<")");

    _startTime = value;
} // setStartTime


// ---------------------------------------------------------------------------------------------------------------------
// Get start time for problem.
double
pylith::problems::TimeDependent::getStartTime(void) const {
    return _startTime;
} // getStartTime


// ---------------------------------------------------------------------------------------------------------------------
// Set end time for problem.
void
pylith::problems::TimeDependent::setEndTime(const double value) {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("endTime(value="<<value<<")");

    if (value < 0.0) {
        std::ostringstream msg;
        msg << "End time (seconds) for problem (" << value << ") must be positive.";
        throw std::runtime_error(msg.str());
    } // if
    _endTime = value;

    PYLITH_METHOD_END;
} // setEndTime


// ---------------------------------------------------------------------------------------------------------------------
// Get end time for problem.
double
pylith::problems::TimeDependent::getEndTime(void) const {
    return _endTime;
} // getEndTime


// ---------------------------------------------------------------------------------------------------------------------
// Set maximum number of time steps.
void
pylith::problems::TimeDependent::setMaxTimeSteps(const size_t value) {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("maxTimeSteps(value="<<value<<")");

    if (value <= 0) {
        std::ostringstream msg;
        msg << "Maximum number of time teps for problem (" << value << ") must be positive.";
        throw std::runtime_error(msg.str());
    } // if
    _maxTimeSteps = value;

    PYLITH_METHOD_END;
} // setMaxTimeSteps


// ---------------------------------------------------------------------------------------------------------------------
// Get maximum number of time steps.
size_t
pylith::problems::TimeDependent::getMaxTimeSteps(void) const {
    return _maxTimeSteps;
} // getMaxTimeSteps


// ---------------------------------------------------------------------------------------------------------------------
// Set initial time step for problem.
void
pylith::problems::TimeDependent::setInitialTimeStep(const double value) {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("setInitialTimeStep(value="<<value<<")");

    if (value < 0.0) {
        std::ostringstream msg;
        msg << "Initial time step (nondimensional) for problem (" << value << ") must be positive.";
        throw std::runtime_error(msg.str());
    } // if
    _dtInitial = value;

    PYLITH_METHOD_END;
} // setInitialTimeStep


// ---------------------------------------------------------------------------------------------------------------------
// Get initial time step for problem.
PetscReal
pylith::problems::TimeDependent::getInitialTimeStep(void) const {
    return _dtInitial;
} // getInitialTimeStep


// ---------------------------------------------------------------------------------------------------------------------
// Set initial conditions.
void
pylith::problems::TimeDependent::setInitialCondition(pylith::problems::InitialCondition* ic[],
                                                     const int numIC) {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("setInitialCondition(ic="<<ic<<", numIC="<<numIC<<")");

    assert( (!ic && 0 == numIC) || (ic && 0 < numIC) );

    _ic.resize(numIC);
    for (int i = 0; i < numIC; ++i) {
        _ic[i] = ic[i];
    } // for

    PYLITH_METHOD_END;
} // setInitialCondition


// ---------------------------------------------------------------------------------------------------------------------
// Should notify observers of solution with initial conditions.
void
pylith::problems::TimeDependent::setShouldNotifyIC(const bool value) {
    _shouldNotifyIC = value;
} // setShouldNotifyIC


// ---------------------------------------------------------------------------------------------------------------------
// Set progress monitor.
void
pylith::problems::TimeDependent::setProgressMonitor(pylith::problems::ProgressMonitorTime* monitor) {
    _monitor = monitor; // :KLUDGE: :TODO: Use shared pointer.
} // setProgressMonitor


// ---------------------------------------------------------------------------------------------------------------------
// Get Petsc DM associated with problem.
PetscDM
pylith::problems::TimeDependent::getPetscDM(void) {
    PYLITH_METHOD_BEGIN;

    PetscDM dm = NULL;
    PetscErrorCode err = TSGetDM(_ts, &dm);PYLITH_CHECK_ERROR(err);

    PYLITH_METHOD_RETURN(dm);
} // getPetscTS


// ---------------------------------------------------------------------------------------------------------------------
/** Get nonlinear solver for problem.
 *
 * @returns PETSc SNES for problem.
 */
PetscSNES
pylith::problems::TimeDependent::getPetscSNES(void) {
    PYLITH_METHOD_BEGIN;

    PetscSNES snes = NULL;
    PetscErrorCode err = TSGetSNES(_ts, &snes);PYLITH_CHECK_ERROR(err);

    PYLITH_METHOD_RETURN(snes);
} // getPetscSNES


// ---------------------------------------------------------------------------------------------------------------------
// Verify configuration.
void
pylith::problems::TimeDependent::verifyConfiguration(void) const {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("Problem::verifyConfiguration(void)");

    Problem::verifyConfiguration();

    assert(_solution);

    // Check to make sure initial conditions are compatible with the solution.
    const size_t numIC = _ic.size();
    for (size_t i = 0; i < numIC; ++i) {
        assert(_ic[i]);
        _ic[i]->verifyConfiguration(*_solution);
    } // for

    PYLITH_METHOD_END;
} // verifyConfiguration


// ---------------------------------------------------------------------------------------------------------------------
// Initialize.
void
pylith::problems::TimeDependent::initialize(void) {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("initialize()");

    Problem::initialize();

    assert(_solution);

    PetscErrorCode err = TSDestroy(&_ts);PYLITH_CHECK_ERROR(err);assert(!_ts);
    const pylith::topology::Mesh& mesh = _solution->mesh();
    err = TSCreate(mesh.comm(), &_ts);PYLITH_CHECK_ERROR(err);assert(_ts);
    err = TSSetType(_ts, TSBEULER);PYLITH_CHECK_ERROR(err); // Backward Euler is default time stepping method.
    // err = TSSetEquationType(_ts, TS_EQ_EXPLICIT);PYLITH_CHECK_ERROR(err);
    err = TSSetExactFinalTime(_ts, TS_EXACTFINALTIME_STEPOVER);PYLITH_CHECK_ERROR(err); // Ok to step over final time.
    err = TSSetFromOptions(_ts);PYLITH_CHECK_ERROR(err);
    err = TSSetApplicationContext(_ts, (void*)this);PYLITH_CHECK_ERROR(err);

    // Set time stepping paramters.
    switch (getSolverType()) {
    case LINEAR:
        PYLITH_COMPONENT_DEBUG("Setting PetscTS problem type to 'linear'.");
        err = TSSetProblemType(_ts, TS_LINEAR);PYLITH_CHECK_ERROR(err);
        break;
    case NONLINEAR:
        PYLITH_COMPONENT_DEBUG("Setting PetscTS problem type to 'nonlinear'.");
        err = TSSetProblemType(_ts, TS_NONLINEAR);PYLITH_CHECK_ERROR(err);
        break;
    default:
        PYLITH_COMPONENT_ERROR("Unknown problem type.");
        throw std::logic_error("Unknown problem type.");
    } // switch
    PYLITH_COMPONENT_DEBUG("Setting PetscTS parameters: "
                           <<"dtInitial="<<_dtInitial
                           <<", startTime="<<_startTime
                           <<", maxTimeSteps="<<_maxTimeSteps
                           <<", endTime="<<_endTime);

    assert(_normalizer);
    const PylithReal timeScale = _normalizer->getTimeScale();
    err = TSSetTime(_ts, _startTime / timeScale);PYLITH_CHECK_ERROR(err);
    err = TSSetTimeStep(_ts, _dtInitial / timeScale);PYLITH_CHECK_ERROR(err);
    err = TSSetMaxSteps(_ts, _maxTimeSteps);PYLITH_CHECK_ERROR(err);
    err = TSSetMaxTime(_ts, _endTime / timeScale);PYLITH_CHECK_ERROR(err);
    err = TSSetDM(_ts, _solution->dmMesh());PYLITH_CHECK_ERROR(err);

    // Set initial solution.
    PYLITH_COMPONENT_DEBUG("Setting PetscTS initial conditions using global vector for solution.");
    _solution->zeroLocal();
    const size_t numIC = _ic.size();
    for (size_t i = 0; i < numIC; ++i) {
        assert(_ic[i]);
        _ic[i]->setValues(_solution, *_normalizer);
    } // for
    _solution->scatterLocalToContext("global");
    err = TSSetSolution(_ts, _solution->scatterVector("global"));PYLITH_CHECK_ERROR(err);
    assert(_observers);
    _observers->setTimeScale(timeScale);

    // Set callbacks.
    PYLITH_COMPONENT_DEBUG("Setting PetscTS callbacks prestep(), poststep(), and computeRHSFunction().");
    err = TSSetPreStep(_ts, prestep);PYLITH_CHECK_ERROR(err);
    err = TSSetPostStep(_ts, poststep);PYLITH_CHECK_ERROR(err);
    err = TSSetRHSFunction(_ts, NULL, computeRHSResidual, (void*)this);PYLITH_CHECK_ERROR(err);

    switch (_formulation) {
    case pylith::problems::Physics::QUASISTATIC:
        PYLITH_COMPONENT_DEBUG("Setting PetscTS callbacks computeRHSJacobian().");
        err = TSSetRHSJacobian(_ts, NULL, NULL, computeRHSJacobian, (void*)this);PYLITH_CHECK_ERROR(err);
    case pylith::problems::Physics::DYNAMIC_IMEX:
        PYLITH_COMPONENT_DEBUG("Setting PetscTS callbacks computeLHSJacobian() and computeLHSFunction().");
        err = TSSetIFunction(_ts, NULL, computeLHSResidual, (void*)this);PYLITH_CHECK_ERROR(err);
        err = TSSetIJacobian(_ts, NULL, NULL, computeLHSJacobian, (void*)this);PYLITH_CHECK_ERROR(err);
    case pylith::problems::Physics::DYNAMIC:
        PYLITH_COMPONENT_DEBUG("Setting up field for inverse of lumped LHS Jacobian.");
        delete _jacobianLHSLumpedInv;_jacobianLHSLumpedInv = new pylith::topology::Field(_solution->mesh());assert(_jacobianLHSLumpedInv);
        _jacobianLHSLumpedInv->cloneSection(*_solution);
        break;
    } // switch

#if 0
    // Set solve type for solution fields defined over the domain (not Lagrange multipliers).
    PetscDS prob = NULL;
    err = DMGetDS(_solution->dmMesh(), &prob);PYLITH_CHECK_ERROR(err);
    PetscInt numFields = 0;
    err = PetscDSGetNumFields(prob, &numFields);PYLITH_CHECK_ERROR(err);
    for (PetscInt iField = 0; iField < numFields; ++iField) {
        err = PetscDSSetImplicit(prob, iField, (_formulationType == IMPLICIT) ? PETSC_TRUE : PETSC_FALSE);
    } // for
#endif
    journal::debug_t debug(pylith::utils::PyreComponent::getName());
    if (debug.state()) {
        PetscDS prob = NULL;
        err = DMGetDS(_solution->dmMesh(), &prob);PYLITH_CHECK_ERROR(err);
        debug << journal::at(__HERE__)
              << "Solution Discretization" << journal::endl;
        PetscDSView(prob, PETSC_VIEWER_STDOUT_SELF);
    } // if

    if (_shouldNotifyIC) {
        _notifyObserversInitialSoln();
    } // if

    if (_monitor) {
        _monitor->open();
    } // if

    PYLITH_METHOD_END;
} // initialize


// ---------------------------------------------------------------------------------------------------------------------
// Solve time-dependent problem.
void
pylith::problems::TimeDependent::solve(void) {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("solve()");

    PetscErrorCode err = TSSolve(_ts, NULL);PYLITH_CHECK_ERROR(err);

    PYLITH_METHOD_END;
} // solve


// ---------------------------------------------------------------------------------------------------------------------
// Perform operations before advancing solution one time step.
void
pylith::problems::TimeDependent::prestep(void) {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("prestep()");

    // Get time of last time step and current time step.
    PetscErrorCode err;
    PylithReal dt;
    PylithReal t;
    err = TSGetTimeStep(_ts, &dt);PYLITH_CHECK_ERROR(err);
    err = TSGetTime(_ts, &t);PYLITH_CHECK_ERROR(err);

    // Prepare constraints for a new time step.
    const size_t numConstraints = _constraints.size();
    for (size_t i = 0; i < numConstraints; ++i) {
        _constraints[i]->prestep(t, dt);
    } // for

    // Prepare integrators for a new time step.
    const size_t numIntegrators = _integrators.size();
    for (size_t i = 0; i < numIntegrators; ++i) {
        _integrators[i]->prestep(t, dt);
    } // for

    PYLITH_METHOD_END;
} // prestep


// ---------------------------------------------------------------------------------------------------------------------
// Perform operations after advancing solution one time step.
void
pylith::problems::TimeDependent::poststep(void) {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("poststep()");

    // Get current solution. After first time step, t==dt, and tindex==1.
    PetscErrorCode err;
    PylithReal t = 0.0, dt = 0.0;
    PylithInt tindex = 0;
    PetscVec solutionVec = NULL;
    err = TSGetTime(_ts, &t);PYLITH_CHECK_ERROR(err);
    err = TSGetTimeStep(_ts, &dt);PYLITH_CHECK_ERROR(err);
    err = TSGetStepNumber(_ts, &tindex);PYLITH_CHECK_ERROR(err);
    err = TSGetSolution(_ts, &solutionVec);PYLITH_CHECK_ERROR(err);

    // Update PyLith view of the solution.
    assert(_solution);
    _solution->scatterVectorToLocal(solutionVec);

    // Update integrators.
    const size_t numIntegrators = _integrators.size();
    for (size_t i = 0; i < numIntegrators; ++i) {
        _integrators[i]->poststep(t, tindex, dt, *_solution);
    } // for

    // Update constraints.
    const size_t numConstraints = _constraints.size();
    for (size_t i = 0; i < numConstraints; ++i) {
        _constraints[i]->poststep(t, tindex, dt, *_solution);
    } // for

    // Notify problem observers of updated solution.
    assert(_observers);
    _observers->notifyObservers(t, tindex, *_solution);

    if (_monitor) {
        assert(_normalizer);
        const PylithReal timeScale = _normalizer->getTimeScale();
        _monitor->update(t*timeScale, _startTime, _endTime);
    } // if

    PYLITH_METHOD_END;
} // poststep


#include <iostream>
// ---------------------------------------------------------------------------------------------------------------------
// Callback static method for computing residual for RHS, G(t,s).
PetscErrorCode
pylith::problems::TimeDependent::computeRHSResidual(PetscTS ts,
                                                    PetscReal t,
                                                    PetscVec solutionVec,
                                                    PetscVec residualVec,
                                                    void* context) {
    PYLITH_METHOD_BEGIN;
    journal::debug_t debug(_TimeDependent::pyreComponent);
    debug << journal::at(__HERE__)
          << "computeRHSResidual(ts="<<ts<<", t="<<t<<", solutionVec="<<solutionVec<<", residualVec="<<residualVec<<", context="<<context<<")" << journal::endl;

    // Get current time step.
    PylithReal dt;
    PetscErrorCode err = TSGetTimeStep(ts, &dt);PYLITH_CHECK_ERROR(err);

    pylith::problems::TimeDependent* problem = (pylith::problems::TimeDependent*)context;
    problem->Problem::computeRHSResidual(residualVec, t, dt, solutionVec);

    if (pylith::problems::Physics::QUASISTATIC != problem->getFormulation()) {
        // Multiply RHS, G(t,s), by M^{-1}
        const PylithReal s_tshift = 1.0; // Keep shift terms on LHS, so use 1.0 for terms moved to RHS.
        problem->Problem::computeLHSJacobianLumpedInv(t, dt, s_tshift, solutionVec);

        assert(problem->_jacobianLHSLumpedInv);
        err = VecPointwiseMult(residualVec, problem->_jacobianLHSLumpedInv->scatterVector("global"), residualVec);PYLITH_CHECK_ERROR(err);
    } // if

    PYLITH_METHOD_RETURN(0);
} // computeRHSResidual


// ---------------------------------------------------------------------------------------------------------------------
// Callback static method for computing Jacobian for RHS, Jacobian of G(t,s).
PetscErrorCode
pylith::problems::TimeDependent::computeRHSJacobian(PetscTS ts,
                                                    PetscReal t,
                                                    PetscVec solutionVec,
                                                    PetscMat jacobianMat,
                                                    PetscMat precondMat,
                                                    void* context) {
    PYLITH_METHOD_BEGIN;
    journal::debug_t debug(_TimeDependent::pyreComponent);
    debug << journal::at(__HERE__)
          << "computeRHSJacobian(ts="<<ts<<", t="<<t<<", solutionVec="<<solutionVec<<", jacobianMat="<<jacobianMat<<", precondMat="<<precondMat<<", context="<<context<<")" << journal::endl;

    // Get current time step.
    PylithReal dt;
    PetscErrorCode err = TSGetTimeStep(ts, &dt);PYLITH_CHECK_ERROR(err);

    pylith::problems::TimeDependent* problem = (pylith::problems::TimeDependent*)context;
    problem->Problem::computeRHSJacobian(jacobianMat, precondMat, t, dt, solutionVec);

    PYLITH_METHOD_RETURN(0);
} // computeRHSJacobian


// ---------------------------------------------------------------------------------------------------------------------
// Callback static method for computing residual for LHS, F(t,s,\dot{s}).
PetscErrorCode
pylith::problems::TimeDependent::computeLHSResidual(PetscTS ts,
                                                    PetscReal t,
                                                    PetscVec solutionVec,
                                                    PetscVec solutionDotVec,
                                                    PetscVec residualVec,
                                                    void* context) {
    PYLITH_METHOD_BEGIN;
    journal::debug_t debug(_TimeDependent::pyreComponent);
    debug << journal::at(__HERE__)
          << "computeLHSResidual(ts="<<ts<<", t="<<t<<", solutionVec="<<solutionVec<<", solutionDotVec="<<solutionDotVec<<", residualVec="<<residualVec<<", context="<<context<<")" << journal::endl;

    // Get current time step.
    PylithReal dt;
    PetscErrorCode err = TSGetTimeStep(ts, &dt);PYLITH_CHECK_ERROR(err);

    pylith::problems::TimeDependent* problem = (pylith::problems::TimeDependent*)context;
    problem->Problem::computeLHSResidual(residualVec, t, dt, solutionVec, solutionDotVec);

    PYLITH_METHOD_RETURN(0);
} // computeLHSResidual


// ---------------------------------------------------------------------------------------------------------------------
// Callback static method for computing Jacobian for LHS, Jacobian of F(t,s,\dot{s}).
PetscErrorCode
pylith::problems::TimeDependent::computeLHSJacobian(PetscTS ts,
                                                    PetscReal t,
                                                    PetscVec solutionVec,
                                                    PetscVec solutionDotVec,
                                                    PetscReal s_tshift,
                                                    PetscMat jacobianMat,
                                                    PetscMat precondMat,
                                                    void* context) {
    PYLITH_METHOD_BEGIN;
    journal::debug_t debug(_TimeDependent::pyreComponent);
    debug << journal::at(__HERE__)
          << "computeLHSJacobian(ts="<<ts<<", t="<<t<<", solutionVec="<<solutionVec<<", solutionDotVec="<<solutionDotVec<<", s_tshift="<<s_tshift<<", jacobianMat="<<jacobianMat<<", precondMat="<<precondMat<<", context="<<context<<")" <<
        journal::endl;

    // Get current time step.
    PylithReal dt;
    PetscErrorCode err = TSGetTimeStep(ts, &dt);PYLITH_CHECK_ERROR(err);

    pylith::problems::TimeDependent* problem = (pylith::problems::TimeDependent*)context;
    problem->Problem::computeLHSJacobian(jacobianMat, precondMat, t, dt, s_tshift, solutionVec, solutionDotVec);

    PYLITH_METHOD_RETURN(0);
} // computeLHSJacobian


// ---------------------------------------------------------------------------------------------------------------------
// Callback static method for operations before advancing solution one time step.
PetscErrorCode
pylith::problems::TimeDependent::prestep(PetscTS ts) {
    PYLITH_METHOD_BEGIN;
    journal::debug_t debug(_TimeDependent::pyreComponent);
    debug << journal::at(__HERE__)
          << "prestep(ts="<<ts<<")" << journal::endl;

    TimeDependent* problem = NULL;
    PetscErrorCode err = TSGetApplicationContext(ts, (void*)&problem);PYLITH_CHECK_ERROR(err);assert(problem);
    problem->prestep();

    PYLITH_METHOD_RETURN(0);
} // prestep


// ---------------------------------------------------------------------------------------------------------------------
// Callback static method for operations after advancing solution one time step.
PetscErrorCode
pylith::problems::TimeDependent::poststep(PetscTS ts) {
    PYLITH_METHOD_BEGIN;
    journal::debug_t debug(_TimeDependent::pyreComponent);
    debug << journal::at(__HERE__)
          << "poststep(ts="<<ts<<")" << journal::endl;

    TimeDependent* problem = NULL;
    PetscErrorCode err = TSGetApplicationContext(ts, (void*)&problem);PYLITH_CHECK_ERROR(err);assert(problem);
    problem->poststep();

    PYLITH_METHOD_RETURN(0);
} // poststep


// ---------------------------------------------------------------------------------------------------------------------
// Notify observers with solution corresponding to initial conditions.
void
pylith::problems::TimeDependent::_notifyObserversInitialSoln(void) {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("_notifyObserversInitialSoln()");

    assert(_normalizer);
    const PylithReal timeScale = _normalizer->getTimeScale();
    const PylithReal tStartNondim = _startTime / timeScale;
    const PylithInt tindex = 0;
    _observers->notifyObservers(tStartNondim, tindex, *_solution);

    const size_t numIntegrators = _integrators.size();
    for (size_t i = 0; i < numIntegrators; ++i) {
        assert(_integrators[i]);
        _integrators[i]->notifyObservers(tStartNondim, tindex, *_solution);
    } // for

    const size_t numConstraints = _constraints.size();
    for (size_t i = 0; i < numConstraints; ++i) {
        assert(_constraints[i]);
        _constraints[i]->notifyObservers(tStartNondim, tindex, *_solution);
    } // for

    PYLITH_METHOD_END;
} // _notifyObserversInitialSoln


// End of file
