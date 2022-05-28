// -*- C++ -*-
//
// ======================================================================
//
// Brad T. Aagaard, U.S. Geological Survey
// Charles A. Williams, GNS Science
// Matthew G. Knepley, University at Buffalo
//
// This code was developed as part of the Computational Infrastructure
// for Geodynamics (http://geodynamics.org).
//
// Copyright (c) 2010-2021 University of California, Davis
//
// See LICENSE.md for license information.
//
// ======================================================================
//

#include <portinfo>

#include "TimeDependent.hh" // implementation of class methods

#include "pylith/feassemble/IntegrationData.hh" // HOLDSA IntegrationData
#include "pylith/topology/Mesh.hh" // USES Mesh
#include "pylith/topology/Field.hh" // USES Field

#include "pylith/feassemble/Integrator.hh" // USES Integrator
#include "pylith/feassemble/Constraint.hh" // USES Constraint
#include "pylith/problems/ObserversSoln.hh" // USES ObserversSoln
#include "pylith/problems/InitialCondition.hh" // USES InitialCondition
#include "pylith/problems/ProgressMonitorTime.hh" // USES ProgressMonitorTime
#include "pylith/utils/PetscOptions.hh" // USES SolverDefaults

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
    _needNewLHSJacobian(true),
    _haveNewLHSJacobian(false),
    _shouldNotifyIC(false) {
    PyreComponent::setName(_TimeDependent::pyreComponent);

    _integrationData->setScalar("dt_jacobian", -1.0);
    _integrationData->setScalar("dt_lhs_jacobian", -1.0);
    _integrationData->setScalar(pylith::feassemble::IntegrationData::t_state, -HUGE_VAL);
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

    _monitor = NULL; // Memory handle in Python. :TODO: Use shared pointer.

    PetscErrorCode err = TSDestroy(&_ts);PYLITH_CHECK_ERROR(err);

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
// Get PETSc time stepper.
PetscTS
pylith::problems::TimeDependent::getPetscTS(void) {
    return _ts;
} // getPetscTS


// ---------------------------------------------------------------------------------------------------------------------
// Verify configuration.
void
pylith::problems::TimeDependent::verifyConfiguration(void) const {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("Problem::verifyConfiguration(void)");

    Problem::verifyConfiguration();

    assert(_integrationData);
    const pylith::topology::Field* solution = _integrationData->getField(pylith::feassemble::IntegrationData::solution);
    assert(solution);

    // Check to make sure initial conditions are compatible with the solution.
    const size_t numIC = _ic.size();
    for (size_t i = 0; i < numIC; ++i) {
        assert(_ic[i]);
        _ic[i]->verifyConfiguration(*solution);
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

    assert(_integrationData);
    pylith::topology::Field* solution = _integrationData->getField(pylith::feassemble::IntegrationData::solution);
    assert(solution);

    PetscErrorCode err = TSDestroy(&_ts);PYLITH_CHECK_ERROR(err);assert(!_ts);
    const pylith::topology::Mesh& mesh = solution->getMesh();
    err = TSCreate(mesh.getComm(), &_ts);PYLITH_CHECK_ERROR(err);assert(_ts);
    err = TSSetType(_ts, TSBEULER);PYLITH_CHECK_ERROR(err); // Backward Euler is default time stepping method.
    err = TSSetExactFinalTime(_ts, TS_EXACTFINALTIME_STEPOVER);PYLITH_CHECK_ERROR(err); // Ok to step over final time.
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
    err = TSSetDM(_ts, solution->getDM());PYLITH_CHECK_ERROR(err);

    // Set initial solution.
    PYLITH_COMPONENT_DEBUG("Setting PetscTS initial conditions using global vector for solution.");
    solution->zeroLocal();
    const size_t numIC = _ic.size();
    for (size_t i = 0; i < numIC; ++i) {
        assert(_ic[i]);
        _ic[i]->setValues(solution, *_normalizer);
    } // for
    PetscVec solutionVector = solution->getGlobalVector();
    solution->scatterLocalToVector(solutionVector);
    err = TSSetSolution(_ts, solutionVector);PYLITH_CHECK_ERROR(err);
    assert(_observers);
    _observers->setTimeScale(timeScale);

    // Initialize residual.
    pylith::topology::Field* residual = new pylith::topology::Field(*solution);assert(residual);
    residual->setLabel("residual");
    _integrationData->setField("residual", residual);

    // Set callbacks.
    PYLITH_COMPONENT_DEBUG("Setting PetscTS callback for poststep().");
    err = TSSetPostStep(_ts, poststep);PYLITH_CHECK_ERROR(err);

    switch (_formulation) {
    case pylith::problems::Physics::QUASISTATIC:
        PYLITH_COMPONENT_DEBUG("Setting PetscTS callbacks computeIFunction() and computeIJacobian().");
        err = TSSetIFunction(_ts, NULL, computeLHSResidual, (void*)this);PYLITH_CHECK_ERROR(err);
        err = TSSetIJacobian(_ts, NULL, NULL, computeLHSJacobian, (void*)this);PYLITH_CHECK_ERROR(err);
        break;
    case pylith::problems::Physics::DYNAMIC_IMEX:
        PYLITH_COMPONENT_DEBUG("Setting PetscTS callbacks computeLHSJacobian() and computeLHSFunction().");
        err = TSSetIFunction(_ts, NULL, computeLHSResidual, (void*)this);PYLITH_CHECK_ERROR(err);
        err = TSSetIJacobian(_ts, NULL, NULL, computeLHSJacobian, (void*)this);PYLITH_CHECK_ERROR(err);
        err = TSSetEquationType(_ts, TS_EQ_EXPLICIT);PYLITH_CHECK_ERROR(err);
    case pylith::problems::Physics::DYNAMIC: {
        PYLITH_COMPONENT_DEBUG("Setting PetscTS callback for computeRHSFunction().");
        err = TSSetRHSFunction(_ts, NULL, computeRHSResidual, (void*)this);PYLITH_CHECK_ERROR(err);
        PYLITH_COMPONENT_DEBUG("Setting up field for inverse of lumped LHS Jacobian.");
        pylith::topology::Field* jacobianLHSLumpedInv = new pylith::topology::Field(*solution);assert(jacobianLHSLumpedInv);
        jacobianLHSLumpedInv->setLabel("JacobianLHS_lumped_inverse");
        jacobianLHSLumpedInv->createGlobalVector();
        _integrationData->setField("jacobian_lhs_lumped_inverse", jacobianLHSLumpedInv);
        break;
    }
    default: {
        PYLITH_COMPONENT_LOGICERROR("Unknown time stepping formulation '" << _formulation << "'.");
    } // default
    } // switch

    pylith::utils::PetscDefaults::set(*solution, _materials[0], _petscDefaults);
    err = TSSetFromOptions(_ts);PYLITH_CHECK_ERROR(err);
    err = TSSetUp(_ts);PYLITH_CHECK_ERROR(err);

#if 0
    // Set solve type for solution fields defined over the domain (not Lagrange multipliers).
    PetscDS prob = NULL;
    err = DMGetDS(solution->getDM(), &prob);PYLITH_CHECK_ERROR(err);
    PetscInt numFields = 0;
    err = PetscDSGetNumFields(prob, &numFields);PYLITH_CHECK_ERROR(err);
    for (PetscInt iField = 0; iField < numFields; ++iField) {
        err = PetscDSSetImplicit(prob, iField, (_formulationType == IMPLICIT) ? PETSC_TRUE : PETSC_FALSE);
    } // for
#endif
    pythia::journal::debug_t debug(pylith::utils::PyreComponent::getName());
    if (debug.state()) {
        PetscDS prob = NULL;
        err = DMGetDS(solution->getDM(), &prob);PYLITH_CHECK_ERROR(err);
        debug << pythia::journal::at(__HERE__)
              << "Solution Discretization" << pythia::journal::endl;
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
    assert(_integrationData);
    pylith::topology::Field* solution = _integrationData->getField(pylith::feassemble::IntegrationData::solution);assert(solution);
    solution->scatterVectorToLocal(solutionVec);
    solution->scatterLocalToOutput();

    // Update integrators.
    const size_t numIntegrators = _integrators.size();
    for (size_t i = 0; i < numIntegrators; ++i) {
        _integrators[i]->poststep(t, tindex, dt, *solution);
    } // for

    // Update constraints.
    const size_t numConstraints = _constraints.size();
    for (size_t i = 0; i < numConstraints; ++i) {
        _constraints[i]->poststep(t, tindex, dt, *solution);
    } // for

    // Notify problem observers of updated solution.
    assert(_observers);
    _observers->notifyObservers(t, tindex, *solution);

    if (_monitor) {
        assert(_normalizer);
        const PylithReal timeScale = _normalizer->getTimeScale();
        _monitor->update(t*timeScale, _startTime, _endTime);
    } // if

    PYLITH_METHOD_END;
} // poststep


// ----------------------------------------------------------------------
// Set solution values according to constraints (Dirichlet BC).
void
pylith::problems::TimeDependent::setSolutionLocal(const PylithReal t,
                                                  PetscVec solutionVec,
                                                  PetscVec solutionDotVec) {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("setSolutionLocal(t="<<t<<", solutionVec="<<solutionVec<<", solutionDotVec="<<solutionDotVec<<")");

    // Update PyLith view of the solution.
    assert(_integrationData);
    pylith::topology::Field* solution = _integrationData->getField(pylith::feassemble::IntegrationData::solution);assert(solution);
    solution->scatterVectorToLocal(solutionVec);
    _integrationData->setScalar(pylith::feassemble::IntegrationData::time, t);

    if (solutionDotVec) {
        pylith::topology::Field* solutionDot = NULL;
        if (!_integrationData->hasField(pylith::feassemble::IntegrationData::solution_dot)) {
            solutionDot = new pylith::topology::Field(*solution);assert(solutionDot);
            solutionDot->setLabel("solutionDot");
            _integrationData->setField(pylith::feassemble::IntegrationData::solution_dot, solutionDot);
        } else {
            solutionDot = _integrationData->getField(pylith::feassemble::IntegrationData::solution_dot);
        } // if/else
        solutionDot->scatterVectorToLocal(solutionDotVec);
    } // if

    const size_t numConstraints = _constraints.size();
    for (size_t i = 0; i < numConstraints; ++i) {
        _constraints[i]->setSolution(_integrationData);
#if 0
        if (solutionDot) {
            _constraints[i]->setSolutionDot(solutionDot, t);
        } // if
#endif
    } // for

    // solution->view("SOLUTION AFTER SETTING VALUES");

    PYLITH_METHOD_END;
} // setSolutionLocal


// ----------------------------------------------------------------------
// Compute RHS residual for G(t,s).
void
pylith::problems::TimeDependent::computeRHSResidual(PetscVec residualVec,
                                                    const PylithReal t,
                                                    const PylithReal dt,
                                                    PetscVec solutionVec) {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("computeRHSResidual(t="<<t<<", dt="<<dt<<", solutionVec="<<solutionVec<<", residualVec="<<residualVec<<")");

    assert(residualVec);
    assert(solutionVec);
    assert(_integrationData);

    if (t != _integrationData->getScalar(pylith::feassemble::IntegrationData::t_state)) { _setState(t); }

    // Update PyLith view of the solution.
    const PetscVec solutionDotVec = NULL;
    setSolutionLocal(t, solutionVec, solutionDotVec);
    _integrationData->setScalar(pylith::feassemble::IntegrationData::time, t);
    _integrationData->setScalar(pylith::feassemble::IntegrationData::time_step, dt);

    // Sum residual contributions across integrators.
    pylith::topology::Field* residual = _integrationData->getField(pylith::feassemble::IntegrationData::residual);assert(residual);
    residual->zeroLocal();
    const size_t numIntegrators = _integrators.size();
    assert(numIntegrators > 0); // must have at least 1 integrator
    for (size_t i = 0; i < numIntegrators; ++i) {
        _integrators[i]->computeRHSResidual(residual, *_integrationData);
    } // for

    // Assemble residual values across processes.
    PetscErrorCode err = VecSet(residualVec, 0.0);PYLITH_CHECK_ERROR(err);
    residual->scatterLocalToVector(residualVec, ADD_VALUES);

    _integrationData->setScalar(pylith::feassemble::IntegrationData::t_state, t);

    PYLITH_METHOD_END;
} // computeRHSResidual


// ----------------------------------------------------------------------
// Compute LHS residual for F(t,s,\dot{s}).
void
pylith::problems::TimeDependent::computeLHSResidual(PetscVec residualVec,
                                                    const PylithReal t,
                                                    const PylithReal dt,
                                                    PetscVec solutionVec,
                                                    PetscVec solutionDotVec) {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("computeLHSResidual(t="<<t<<", dt="<<dt<<", solutionVec="<<solutionVec<<", solutionDotVec="<<solutionDotVec<<", residualVec="<<residualVec<<")");

    assert(residualVec);
    assert(solutionVec);
    assert(solutionDotVec);
    assert(_integrationData);

    if (t != _integrationData->getScalar(pylith::feassemble::IntegrationData::t_state)) { _setState(t); }

    // Update PyLith view of the solution.
    setSolutionLocal(t, solutionVec, solutionDotVec);
    _integrationData->setScalar(pylith::feassemble::IntegrationData::time, t);
    _integrationData->setScalar(pylith::feassemble::IntegrationData::time_step, dt);

    // Sum residual across integrators.
    pylith::topology::Field* residual = _integrationData->getField(pylith::feassemble::IntegrationData::residual);assert(residual);
    residual->zeroLocal();
    const int numIntegrators = _integrators.size();
    assert(numIntegrators > 0); // must have at least 1 integrator
    for (int i = 0; i < numIntegrators; ++i) {
        _integrators[i]->computeLHSResidual(residual, *_integrationData);
    } // for

    // Assemble residual values across processes.
    PetscErrorCode err = VecSet(residualVec, 0.0);PYLITH_CHECK_ERROR(err);
    residual->scatterLocalToVector(residualVec, ADD_VALUES);

    _integrationData->setScalar(pylith::feassemble::IntegrationData::t_state, t);

    PYLITH_METHOD_END;
} // computeLHSResidual


// ----------------------------------------------------------------------
// Compute LHS Jacobian for F(t,s,\dot{s}) for implicit time stepping.
void
pylith::problems::TimeDependent::computeLHSJacobian(PetscMat jacobianMat,
                                                    PetscMat precondMat,
                                                    const PylithReal t,
                                                    const PylithReal dt,
                                                    const PylithReal s_tshift,
                                                    PetscVec solutionVec,
                                                    PetscVec solutionDotVec) {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("computeLHSJacobian(t="<<t<<", dt="<<dt<<", s_tshift="<<s_tshift<<", solutionVec="<<solutionVec<<", solutionDotVec="<<solutionDotVec<<", jacobianMat="<<jacobianMat<<", precondMat="<<precondMat<<")");

    assert(jacobianMat);
    assert(precondMat);
    assert(solutionVec);
    assert(solutionDotVec);
    assert(s_tshift > 0);

    if (!_needNewJacobian(dt)) {
        PYLITH_COMPONENT_DEBUG("KEEP LHS Jacobian; t=" << t << ", dt=" << dt);
        _haveNewLHSJacobian = false;
        PYLITH_METHOD_END;
    } // if
    PYLITH_COMPONENT_DEBUG("NEW LHS Jacobian; t=" << t << ", dt=" << dt);

    assert(_integrationData);
    const pylith::topology::Field* solution = _integrationData->getField(pylith::feassemble::IntegrationData::solution);assert(solution);
    _integrationData->setScalar(pylith::feassemble::IntegrationData::time, t);
    _integrationData->setScalar(pylith::feassemble::IntegrationData::time_step, dt);
    _integrationData->setScalar(pylith::feassemble::IntegrationData::s_tshift, s_tshift);

    // Zero LHS Jacobian
    PetscErrorCode err = 0;
    PetscDS solnDS = NULL;
    PetscBool hasJacobian = PETSC_FALSE;
    err = DMGetDS(solution->getDM(), &solnDS);PYLITH_CHECK_ERROR(err);
    err = PetscDSHasJacobian(solnDS, &hasJacobian);PYLITH_CHECK_ERROR(err);
    if (hasJacobian) { err = MatZeroEntries(jacobianMat);PYLITH_CHECK_ERROR(err); }
    err = MatZeroEntries(precondMat);PYLITH_CHECK_ERROR(err);

    // Update PyLith view of the solution.
    setSolutionLocal(t, solutionVec, solutionDotVec);

    // Sum Jacobian contributions across integrators.
    const size_t numIntegrators = _integrators.size();
    for (size_t i = 0; i < numIntegrators; ++i) {
        _integrators[i]->computeLHSJacobian(jacobianMat, precondMat, *_integrationData);
    } // for

    _needNewLHSJacobian = false;
    _haveNewLHSJacobian = true;

    _integrationData->setScalar(pylith::feassemble::IntegrationData::dt_jacobian, dt);

    // Solver handles assembly.

    PYLITH_METHOD_END;
} // computeLHSJacobian


// ----------------------------------------------------------------------
// Compute inverse of LHS Jacobian for F(t,s,\dot{s}) for explicit time stepping.
void
pylith::problems::TimeDependent::computeLHSJacobianLumpedInv(const PylithReal t,
                                                             const PylithReal dt,
                                                             const PylithReal s_tshift,
                                                             PetscVec solutionVec) {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("TimeDependent::computeLHSJacobianLumpedInv(t="<<t<<", dt="<<dt<<", s_tshift="<<s_tshift<<", solutionVec="<<solutionVec<<")");

    assert(solutionVec);
    assert(s_tshift > 0);
    assert(_integrationData);

    const size_t numIntegrators = _integrators.size();

    // Check to see if we need to compute LHS Jacobian.
    bool needNewLHSJacobianLumped = false;
    const bool dtChanged = dt != _integrationData->getScalar(pylith::feassemble::IntegrationData::dt_lumped_jacobian_inverse);
    for (size_t i = 0; i < numIntegrators; ++i) {
        if (_integrators[i]->needNewLHSJacobianLumped(dtChanged)) {
            needNewLHSJacobianLumped = true;
            break;
        } // if
    } // for
    if (!needNewLHSJacobianLumped) { PYLITH_METHOD_END; }

    // Set jacobian to zero.
    pylith::topology::Field* jacobianLumpedInv = _integrationData->getField(pylith::feassemble::IntegrationData::lumped_jacobian_inverse);
    jacobianLumpedInv->zeroLocal();

    // Update PyLith view of the solution.
    const PetscVec solutionDotVec = NULL;
    setSolutionLocal(t, solutionVec, solutionDotVec);
    _integrationData->setScalar(pylith::feassemble::IntegrationData::time, t);
    _integrationData->setScalar(pylith::feassemble::IntegrationData::time_step, dt);
    _integrationData->setScalar(pylith::feassemble::IntegrationData::s_tshift, s_tshift);

    // Sum Jacobian contributions across integrators.
    for (size_t i = 0; i < numIntegrators; ++i) {
        _integrators[i]->computeLHSJacobianLumpedInv(jacobianLumpedInv, *_integrationData);
    } // for

    // Insert values into global vector.
    jacobianLumpedInv->scatterLocalToVector(jacobianLumpedInv->getGlobalVector());

    _integrationData->setScalar(pylith::feassemble::IntegrationData::dt_lumped_jacobian_inverse, dt);

    PYLITH_METHOD_END;
} // computeLHSJacobianLumpedInv


// ---------------------------------------------------------------------------------------------------------------------
// Callback static method for computing residual for RHS, G(t,s).
PetscErrorCode
pylith::problems::TimeDependent::computeRHSResidual(PetscTS ts,
                                                    PetscReal t,
                                                    PetscVec solutionVec,
                                                    PetscVec residualVec,
                                                    void* context) {
    PYLITH_METHOD_BEGIN;
    pythia::journal::debug_t debug(_TimeDependent::pyreComponent);
    debug << pythia::journal::at(__HERE__)
          << "computeRHSResidual(ts="<<ts<<", t="<<t<<", solutionVec="<<solutionVec<<", residualVec="<<residualVec<<", context="<<context<<")" << pythia::journal::endl;

    // Get current time step.
    PylithReal dt;
    PetscErrorCode err = TSGetTimeStep(ts, &dt);PYLITH_CHECK_ERROR(err);

    pylith::problems::TimeDependent* problem = (pylith::problems::TimeDependent*)context;assert(problem);
    assert(problem->_integrationData);
    const bool hasLumpedJacobianInverse = problem->_integrationData->hasField(pylith::feassemble::IntegrationData::lumped_jacobian_inverse);
    if (hasLumpedJacobianInverse) {
        const PylithReal s_tshift = 1.0; // Keep shift terms on LHS, so use 1.0 for terms moved to RHS.
        problem->computeLHSJacobianLumpedInv(t, dt, s_tshift, solutionVec);
    } // if

    problem->computeRHSResidual(residualVec, t, dt, solutionVec);

    if (hasLumpedJacobianInverse) {
        // Multiply RHS, G(t,s), by M^{-1}
        const pylith::topology::Field* jacobianLumpedInv =
            problem->_integrationData->getField(pylith::feassemble::IntegrationData::lumped_jacobian_inverse);assert(jacobianLumpedInv);
        err = VecPointwiseMult(residualVec, jacobianLumpedInv->getGlobalVector(), residualVec);PYLITH_CHECK_ERROR(err);
    } // if

    PYLITH_METHOD_RETURN(0);
} // computeRHSResidual


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
    pythia::journal::debug_t debug(_TimeDependent::pyreComponent);
    debug << pythia::journal::at(__HERE__)
          << "computeLHSResidual(ts="<<ts<<", t="<<t<<", solutionVec="<<solutionVec<<", solutionDotVec="<<solutionDotVec<<", residualVec="<<residualVec<<", context="<<context<<")" << pythia::journal::endl;

    // Get current time step.
    PylithReal dt;
    PetscErrorCode err = TSGetTimeStep(ts, &dt);PYLITH_CHECK_ERROR(err);

    pylith::problems::TimeDependent* problem = (pylith::problems::TimeDependent*)context;
    problem->computeLHSResidual(residualVec, t, dt, solutionVec, solutionDotVec);

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
    pythia::journal::debug_t debug(_TimeDependent::pyreComponent);
    debug << pythia::journal::at(__HERE__)
          << "computeLHSJacobian(ts="<<ts<<", t="<<t<<", solutionVec="<<solutionVec<<", solutionDotVec="<<solutionDotVec<<", s_tshift="<<s_tshift<<", jacobianMat="<<jacobianMat<<", precondMat="<<precondMat<<", context="<<context<<")" <<
        pythia::journal::endl;

    // Get current time step.
    PylithReal dt;
    PetscErrorCode err = TSGetTimeStep(ts, &dt);PYLITH_CHECK_ERROR(err);

    pylith::problems::TimeDependent* problem = (pylith::problems::TimeDependent*)context;
    problem->computeLHSJacobian(jacobianMat, precondMat, t, dt, s_tshift, solutionVec, solutionDotVec);

    PYLITH_METHOD_RETURN(0);
} // computeLHSJacobian


// ---------------------------------------------------------------------------------------------------------------------
// Callback static method for operations after advancing solution one time step.
PetscErrorCode
pylith::problems::TimeDependent::poststep(PetscTS ts) {
    PYLITH_METHOD_BEGIN;
    pythia::journal::debug_t debug(_TimeDependent::pyreComponent);
    debug << pythia::journal::at(__HERE__)
          << "poststep(ts="<<ts<<")" << pythia::journal::endl;

    TimeDependent* problem = NULL;
    PetscErrorCode err = TSGetApplicationContext(ts, (void*)&problem);PYLITH_CHECK_ERROR(err);assert(problem);
    problem->poststep();

    PYLITH_METHOD_RETURN(0);
} // poststep


// ---------------------------------------------------------------------------------------------------------------------
// Check whether we need to reform the Jacobian.
bool
pylith::problems::TimeDependent::_needNewJacobian(const PylithReal dt) {
    PYLITH_METHOD_BEGIN;

    // If we already know we need to recompute the LHS Jacobian, then return true.
    if (_needNewLHSJacobian) { PYLITH_METHOD_RETURN(true); }

    assert(_integrationData);
    const bool dtChanged = dt != _integrationData->getScalar(pylith::feassemble::IntegrationData::dt_jacobian);
    const size_t numIntegrators = _integrators.size();

    for (size_t i = 0; i < numIntegrators; ++i) {
        if (_integrators[i]->needNewLHSJacobian(dtChanged)) {
            _needNewLHSJacobian = true;
            break;
        } // if
    } // for

    PYLITH_METHOD_RETURN(_needNewLHSJacobian);
} // _needNewJacobian


// ---------------------------------------------------------------------------------------------------------------------
// Set state (auxiliary field values) of system for time t.
void
pylith::problems::TimeDependent::_setState(const PylithReal t) {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("_setState(t="<<t<<")");

    // Update constraint values to current time, t.
    const size_t numConstraints = _constraints.size();
    for (size_t i = 0; i < numConstraints; ++i) {
        _constraints[i]->setState(t);
    } // for

    // Prepare integrators for a new time step.
    const size_t numIntegrators = _integrators.size();
    for (size_t i = 0; i < numIntegrators; ++i) {
        _integrators[i]->setState(t);
    } // for

    PYLITH_METHOD_END;
} // _setState


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

    assert(_integrationData);
    const pylith::topology::Field* solution = _integrationData->getField(pylith::feassemble::IntegrationData::solution);assert(solution);
    _observers->notifyObservers(tStartNondim, tindex, *solution);

    const size_t numIntegrators = _integrators.size();
    for (size_t i = 0; i < numIntegrators; ++i) {
        assert(_integrators[i]);
        _integrators[i]->notifyObservers(tStartNondim, tindex, *solution);
    } // for

    const size_t numConstraints = _constraints.size();
    for (size_t i = 0; i < numConstraints; ++i) {
        assert(_constraints[i]);
        _constraints[i]->notifyObservers(tStartNondim, tindex, *solution);
    } // for

    PYLITH_METHOD_END;
} // _notifyObserversInitialSoln


// End of file
