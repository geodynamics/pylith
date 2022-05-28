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

#include "GreensFns.hh" // implementation of class methods

#include "pylith/feassemble/IntegrationData.hh" // HOLDSA IntegrationData
#include "pylith/topology/Mesh.hh" // USES Mesh
#include "pylith/topology/Field.hh" // USES Field

#include "pylith/faults/FaultCohesiveImpulses.hh" // USES FaultCohesiveImpulses
#include "pylith/feassemble/IntegratorInterface.hh" // USES IntegratorInterface
#include "pylith/feassemble/Constraint.hh" // USES Constraint
#include "pylith/problems/ObserversSoln.hh" // USES ObserversSoln
#include "pylith/problems/ProgressMonitorStep.hh" // USES ProgressMonitorStep
#include "pylith/utils/PetscOptions.hh" // USES SolverDefaults

#include "spatialdata/units/Nondimensional.hh" // USES Nondimensional

#include "petscsnes.h" // USES PetscSNES

#include "pylith/utils/error.hh" // USES PYLITH_CHECK_ERROR
#include "pylith/utils/journals.hh" // USES PYLITH_COMPONENT_*
#include <cassert> // USES assert()

// ------------------------------------------------------------------------------------------------
namespace pylith {
    namespace problems {
        class _GreensFns {
public:

            static const char* pyreComponent;
        }; // _GreensFns

        const char* _GreensFns::pyreComponent = "greensfns";
    } // problems
} // pylith

// ------------------------------------------------------------------------------------------------
// Constructor
pylith::problems::GreensFns::GreensFns(void) :
    _faultLabelName(pylith::topology::Mesh::cells_label_name),
    _faultLabelValue(100),
    _faultImpulses(NULL),
    _integratorImpulses(NULL),
    _snes(NULL),
    _monitor(NULL) {
    PyreComponent::setName(_GreensFns::pyreComponent);

    _integrationData->setScalar("dt_jacobian", -1.0);
    _integrationData->setScalar("dt_lhs_jacobian", -1.0);
    _integrationData->setScalar(pylith::feassemble::IntegrationData::t_state, -HUGE_VAL);
} // constructor


// ------------------------------------------------------------------------------------------------
// Destructor
pylith::problems::GreensFns::~GreensFns(void) {
    deallocate();
} // destructor


// ------------------------------------------------------------------------------------------------
// Deallocate PETSc and local data structures.
void
pylith::problems::GreensFns::deallocate(void) {
    PYLITH_METHOD_BEGIN;

    Problem::deallocate();

    _faultImpulses = NULL; // Memory handle in Python. :TODO: Use shared pointer.
    _integratorImpulses = NULL; // Memory handle in Problem. :TODO: Use shared pointer.

    _monitor = NULL; // Memory handle in Python. :TODO: Use shared pointer.

    PetscErrorCode err = SNESDestroy(&_snes);PYLITH_CHECK_ERROR(err);

    PYLITH_METHOD_END;
} // deallocate


// ------------------------------------------------------------------------------------------------
// Set name of label for fault with impulses.
void
pylith::problems::GreensFns::setFaultLabelName(const char* value) {
    PYLITH_COMPONENT_DEBUG("setFaultLabelName(value="<<value<<")");

    if (strlen(value) == 0) {
        throw std::runtime_error("Empty string given for name of label for fault with impulses.");
    } // if

    _faultLabelName = value;
} // setFaultLabelName


// ------------------------------------------------------------------------------------------------
// Get label name for fault with impulses.
const char*
pylith::problems::GreensFns::getFaultLabelName(void) const {
    return _faultLabelName.c_str();
} // getFaultLabelName


// ------------------------------------------------------------------------------------------------
// Set value of label for fault with impulses.
void
pylith::problems::GreensFns::setFaultLabelValue(const int value) {
    PYLITH_COMPONENT_DEBUG("fsetFaultLabelValue(value="<<value<<")");

    _faultLabelValue = value;
} // setFaultLabelvalue


// ------------------------------------------------------------------------------------------------
// Get label value for fault with impulses.
int
pylith::problems::GreensFns::getFaultLabelValue(void) const {
    return _faultLabelValue;
} // getFaultLabelValue


// ------------------------------------------------------------------------------------------------
// Set progress monitor.
void
pylith::problems::GreensFns::setProgressMonitor(pylith::problems::ProgressMonitorStep* monitor) {
    _monitor = monitor; // :KLUDGE: :TODO: Use shared pointer.
} // setProgressMonitor


// ------------------------------------------------------------------------------------------------
// Get Petsc DM associated with problem.
PetscDM
pylith::problems::GreensFns::getPetscDM(void) {
    PYLITH_METHOD_BEGIN;

    PetscDM dm = NULL;
    PetscErrorCode err = SNESGetDM(_snes, &dm);PYLITH_CHECK_ERROR(err);

    PYLITH_METHOD_RETURN(dm);
} // getPetscDM


// ------------------------------------------------------------------------------------------------
// Verify configuration.
void
pylith::problems::GreensFns::verifyConfiguration(void) const {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("Problem::verifyConfiguration(void)");

    Problem::verifyConfiguration();

    // Verify we have fault for the impulses.
    const size_t numInterfaces = _interfaces.size();
    pylith::faults::FaultCohesiveImpulses* faultImpulses = NULL;
    for (size_t i = 0; i < numInterfaces; ++i) {
        if ((_faultLabelName == std::string(_interfaces[i]->getSurfaceLabelName())) &&
            (_faultLabelValue == _interfaces[i]->getSurfaceLabelValue())) {
            faultImpulses = dynamic_cast<pylith::faults::FaultCohesiveImpulses*>(_interfaces[i]);
            if (!faultImpulses) {
                std::ostringstream msg;
                msg << "Found fault with "<<_faultLabelName<<"="<<_faultLabelValue
                    <<" in interfaces for imposing impulses, but type is not FaultCohesiveImpulses.";
                throw std::runtime_error(msg.str());
            } // if
            break;
        } // if
    } // for
    if (!faultImpulses) {
        std::ostringstream msg;
        msg << "Could not find fault with "<<_faultLabelName<<"="<<_faultLabelValue<<" in interfaces for imposing impulses.";
        throw std::runtime_error(msg.str());
    } // if

    PYLITH_METHOD_END;
} // verifyConfiguration


// ------------------------------------------------------------------------------------------------
// Initialize.
void
pylith::problems::GreensFns::initialize(void) {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("initialize()");

    Problem::initialize();

    assert(_integrationData);

    // Find fault on which to put the impulses.
    assert(!_faultImpulses);
    const size_t numInterfaces = _interfaces.size();
    for (size_t i = 0; i < numInterfaces; ++i) {
        if ((_faultLabelName == std::string(_interfaces[i]->getSurfaceLabelName())) &&
            (_faultLabelValue == _interfaces[i]->getSurfaceLabelValue())) {
            _faultImpulses = dynamic_cast<pylith::faults::FaultCohesiveImpulses*>(_interfaces[i]);
        } // if
    } // for
    assert(_faultImpulses);

    PetscErrorCode err = SNESDestroy(&_snes);PYLITH_CHECK_ERROR(err);assert(!_snes);
    assert(_integrationData);
    pylith::topology::Field* solution = _integrationData->getField(pylith::feassemble::IntegrationData::solution);
    assert(solution);

    err = SNESCreate(solution->getMesh().getComm(), &_snes);PYLITH_CHECK_ERROR(err);assert(_snes);
    err = SNESSetApplicationContext(_snes, (void*)this);PYLITH_CHECK_ERROR(err);
    err = SNESSetDM(_snes, solution->getDM());PYLITH_CHECK_ERROR(err);

    PetscVec solutionVector = solution->getGlobalVector();
    solution->scatterLocalToVector(solutionVector);
    err = SNESSetSolution(_snes, solutionVector);PYLITH_CHECK_ERROR(err);

    // Initialize solution_dot.
    pylith::topology::Field* solutionDot = new pylith::topology::Field(*solution);assert(solutionDot);
    solutionDot->setLabel("solution_dot");
    _integrationData->setField(pylith::feassemble::IntegrationData::solution_dot, solutionDot);

    // Initialize residual.
    pylith::topology::Field* residual = new pylith::topology::Field(*solution);assert(residual);
    residual->setLabel("residual");
    _integrationData->setField(pylith::feassemble::IntegrationData::residual, residual);

    _integrationData->setScalar(pylith::feassemble::IntegrationData::time, 0.0);
    _integrationData->setScalar(pylith::feassemble::IntegrationData::time_step, 1.0);
    _integrationData->setScalar(pylith::feassemble::IntegrationData::s_tshift, 0.0);

    switch (_formulation) {
    case pylith::problems::Physics::QUASISTATIC:
        PYLITH_COMPONENT_DEBUG("Setting PetscSNES callbacks SNESSetFunction() and SNESSetJacobian().");
        err = SNESSetFunction(_snes, NULL, computeResidual, (void*)this);PYLITH_CHECK_ERROR(err);
        err = SNESSetJacobian(_snes, NULL, NULL, computeJacobian, (void*)this);PYLITH_CHECK_ERROR(err);
        err = SNESSetType(_snes, SNESKSPONLY);PYLITH_CHECK_ERROR(err);
        break;
    case pylith::problems::Physics::DYNAMIC_IMEX:
        PYLITH_COMPONENT_LOGICERROR("Dynamic Green's functions problems not yet supported.");
        break;
    case pylith::problems::Physics::DYNAMIC:
        PYLITH_COMPONENT_LOGICERROR("Dynamic Green's functions problems not yet supported.");
        break;
    default: {
        PYLITH_COMPONENT_LOGICERROR("Unknown Green's functions formulation '" << _formulation << "'.");
    } // default
    } // switch

    pylith::utils::PetscDefaults::set(*solution, _materials[0], _petscDefaults);
    err = SNESSetFromOptions(_snes);PYLITH_CHECK_ERROR(err);
    err = SNESSetUp(_snes);PYLITH_CHECK_ERROR(err);

    // Get integrator for fault with impulses.
    assert(!_integratorImpulses);
    const size_t numIntegrators = _integrators.size();
    for (size_t i = 0; i < numIntegrators; ++i) {
        pylith::feassemble::IntegratorInterface* integrator = dynamic_cast<pylith::feassemble::IntegratorInterface*>(_integrators[i]);
        if (integrator && (_faultLabelName == std::string(integrator->getSurfaceLabelName()))) {
            _integratorImpulses = _integrators[i];
        } // if
    } // for
    if (!_integratorImpulses) {
        std::ostringstream msg;
        msg << "Could not find integrator for fault "<<_faultLabelName<<"="<<_faultLabelValue<<" in integrators for problem.";
        throw std::runtime_error(msg.str());
    } // if

    pythia::journal::debug_t debug(pylith::utils::PyreComponent::getName());
    if (debug.state()) {
        PetscDS dsSoln = NULL;
        err = DMGetDS(solution->getDM(), &dsSoln);PYLITH_CHECK_ERROR(err);
        debug << pythia::journal::at(__HERE__)
              << "Solution Discretization" << pythia::journal::endl;
        PetscDSView(dsSoln, PETSC_VIEWER_STDOUT_SELF);
    } // if

    if (_monitor) {
        _monitor->open();
    } // if

    PYLITH_METHOD_END;
} // initialize


// ------------------------------------------------------------------------------------------------
// Solve Green's functions problem.
void
pylith::problems::GreensFns::solve(void) {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("solve()");

    assert(_integrationData);
    pylith::topology::Field* solution = _integrationData->getField(pylith::feassemble::IntegrationData::solution);
    assert(solution);
    pylith::topology::Field* residual = _integrationData->getField(pylith::feassemble::IntegrationData::residual);assert(residual);

    const size_t numImpulses = _faultImpulses->getNumImpulses();
    PYLITH_COMPONENT_DEBUG("Using " << numImpulses << " impulses for Green's functions.");

    int mpiRank = 0;
    MPI_Comm_rank(PETSC_COMM_WORLD, &mpiRank);
    const PylithReal tolerance = 1.0e-4;
    for (size_t i = 0; i < numImpulses; ++i) {
        if (0 == mpiRank) {
            PYLITH_COMPONENT_INFO("Computing Green's function " << i+1 << " of " << numImpulses << ".");
        } // if

        // Update impulse on fault
        const PetscReal impulseReal = i + tolerance;
        _integratorImpulses->setState(impulseReal);

        PetscErrorCode err = SNESSolve(_snes, residual->getGlobalVector(), solution->getGlobalVector());PYLITH_CHECK_ERROR(err);

        solution->scatterVectorToLocal(solution->getGlobalVector());
        solution->scatterLocalToOutput();
        poststep(size_t(impulseReal));
    } // for

    PYLITH_METHOD_END;
} // solve


// ------------------------------------------------------------------------------------------------
// Perform operations after advancing solution of one impulse.
void
pylith::problems::GreensFns::poststep(const size_t impulse) {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("poststep(impulse"<<impulse<<")");

    // Get current solution. ParaView needs valid times.
    const PetscReal t = impulse / _normalizer->getTimeScale();
    const PetscReal dt = 1.0;

    assert(_integrationData);
    pylith::topology::Field* solution = _integrationData->getField(pylith::feassemble::IntegrationData::solution);
    assert(solution);

    // Update integrators.
    const size_t numIntegrators = _integrators.size();
    for (size_t i = 0; i < numIntegrators; ++i) {
        _integrators[i]->poststep(t, impulse, dt, *solution);
    } // for

    // Update constraints.
    const size_t numConstraints = _constraints.size();
    for (size_t i = 0; i < numConstraints; ++i) {
        _constraints[i]->poststep(t, impulse, dt, *solution);
    } // for

    // Notify problem observers of updated solution.
    assert(_observers);
    _observers->notifyObservers(t, impulse, *solution);

    // Get number of impulses for monitor
    const size_t numImpulses = _faultImpulses->getNumImpulses();
    if (_monitor) {
        assert(_normalizer);
        _monitor->update(impulse, 0, numImpulses);
    } // if

    PYLITH_METHOD_END;
} // poststep


// ------------------------------------------------------------------------------------------------
// Set solution values according to constraints (Dirichlet BC).
void
pylith::problems::GreensFns::setSolutionLocal(PetscVec solutionVec) {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("setSolutionLocal(solutionVec="<<solutionVec<<")");

    // Update PyLith view of the solution.
    assert(_integrationData);
    pylith::topology::Field* solution = _integrationData->getField(pylith::feassemble::IntegrationData::solution);assert(solution);
    solution->scatterVectorToLocal(solutionVec);

    const size_t numConstraints = _constraints.size();
    for (size_t i = 0; i < numConstraints; ++i) {
        _constraints[i]->setSolution(_integrationData);
    } // for

    // solution->view("SOLUTION AFTER SETTING VALUES");

    PYLITH_METHOD_END;
} // setSolutionLocal


// ------------------------------------------------------------------------------------------------
// Compute residual.
void
pylith::problems::GreensFns::computeResidual(PetscVec residualVec,
                                             PetscVec solutionVec) {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("computeResidual(solutionVec="<<solutionVec<<", residualVec="<<residualVec<<")");

    assert(residualVec);
    assert(solutionVec);
    assert(_integrationData);
    pylith::topology::Field* residual = _integrationData->getField(pylith::feassemble::IntegrationData::residual);
    assert(residual);

    // Update PyLith view of the solution.
    setSolutionLocal(solutionVec);

    // Sum residual across integrators.
    residual->zeroLocal();
    const int numIntegrators = _integrators.size();
    assert(numIntegrators > 0); // must have at least 1 integrator
    for (int i = 0; i < numIntegrators; ++i) {
        _integrators[i]->computeLHSResidual(residual, *_integrationData);
    } // for

    // Assemble residual values across processes.
    PetscErrorCode err = VecSet(residualVec, 0.0);PYLITH_CHECK_ERROR(err);
    residual->scatterLocalToVector(residualVec, ADD_VALUES);

    PYLITH_METHOD_END;
} // computeResidual


// ------------------------------------------------------------------------------------------------
// Compute Jacobian.
void
pylith::problems::GreensFns::computeJacobian(PetscMat jacobianMat,
                                             PetscMat precondMat,
                                             PetscVec solutionVec) {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("GreensFns::computeJacobian(solutionVec="<<solutionVec<<",jacobianMat="<<jacobianMat<<", precondMat="<<precondMat<<")");

    assert(jacobianMat);
    assert(solutionVec);

    assert(_integrationData);
    pylith::topology::Field* solution = _integrationData->getField(pylith::feassemble::IntegrationData::solution);
    assert(solution);

    // Zero Jacobian
    PetscErrorCode err = 0;
    PetscDS dsSoln = NULL;
    err = DMGetDS(solution->getDM(), &dsSoln);PYLITH_CHECK_ERROR(err);

    PetscBool hasJacobian = PETSC_FALSE;
    err = PetscDSHasJacobian(dsSoln, &hasJacobian);PYLITH_CHECK_ERROR(err);
    if (hasJacobian) { err = MatZeroEntries(jacobianMat);PYLITH_CHECK_ERROR(err); }

    PetscBool hasPreconditioner = PETSC_FALSE;
    err = PetscDSHasJacobianPreconditioner(dsSoln, &hasPreconditioner);PYLITH_CHECK_ERROR(err);
    if (hasPreconditioner) { err = MatZeroEntries(precondMat);PYLITH_CHECK_ERROR(err); }

    // Update PyLith view of the solution.
    setSolutionLocal(solutionVec);

    // Sum Jacobian contributions across integrators.
    const size_t numIntegrators = _integrators.size();
    for (size_t i = 0; i < numIntegrators; ++i) {
        _integrators[i]->computeLHSJacobian(jacobianMat, precondMat, *_integrationData);
    } // for

    // Solver handles assembly.

    PYLITH_METHOD_END;
} // computeJacobian


// ------------------------------------------------------------------------------------------------
// Callback static method for computing residual.
PetscErrorCode
pylith::problems::GreensFns::computeResidual(PetscSNES snes,
                                             PetscVec solutionVec,
                                             PetscVec residualVec,
                                             void* context) {
    PYLITH_METHOD_BEGIN;
    pythia::journal::debug_t debug(_GreensFns::pyreComponent);
    debug << pythia::journal::at(__HERE__)
          << "computeResidual(snes="<<snes<<", solutionVec="<<solutionVec<<", residualVec="<<residualVec<<", context="<<context<<")" << pythia::journal::endl;

    pylith::problems::GreensFns* problem = (pylith::problems::GreensFns*)context;
    problem->computeResidual(residualVec, solutionVec);

    PYLITH_METHOD_RETURN(0);
} // computeResidual


// ------------------------------------------------------------------------------------------------
// Callback static method for computing Jacobian.
PetscErrorCode
pylith::problems::GreensFns::computeJacobian(PetscSNES snes,
                                             PetscVec solutionVec,
                                             PetscMat jacobianMat,
                                             PetscMat precondMat,
                                             void* context) {
    PYLITH_METHOD_BEGIN;
    pythia::journal::debug_t debug(_GreensFns::pyreComponent);
    debug << pythia::journal::at(__HERE__)
          << "computeJacobian(snes="<<snes<<", solutionVec="<<solutionVec<<", jacobianMat="<<jacobianMat<<", precondMat="<<precondMat<<", context="<<context<<")" <<
        pythia::journal::endl;

    pylith::problems::GreensFns* problem = (pylith::problems::GreensFns*)context;
    problem->computeJacobian(jacobianMat, precondMat, solutionVec);

    PYLITH_METHOD_RETURN(0);
} // computeJacobian


// End of file
