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

#include "pylith/topology/Mesh.hh" // USES Mesh
#include "pylith/topology/Field.hh" // USES Field
#include "pylith/faults/FaultCohesiveImpulses.hh" // USES FaultCohesiveImpulses
#include "pylith/feassemble/Integrator.hh" // USES Integrator
#include "pylith/feassemble/Constraint.hh" // USES Constraint
#include "pylith/problems/ObserversSoln.hh" // USES ObserversSoln
#include "pylith/problems/ProgressMonitorStep.hh" // USES ProgressMonitorStep

#include "spatialdata/units/Nondimensional.hh" // USES Nondimensional

#include "petscsnes.h" // USES PetscSNES

#include "pylith/utils/error.hh" // USES PYLITH_CHECK_ERROR
#include "pylith/utils/journals.hh" // USES PYLITH_COMPONENT_*
#include <cassert> // USES assert()

// ---------------------------------------------------------------------------------------------------------------------
namespace pylith {
    namespace problems {
        class _GreensFns {
public:

            static const char* pyreComponent;
        }; // _GreensFns

        const char* _GreensFns::pyreComponent = "GreensFns";
    } // problems
} // pylith

// ---------------------------------------------------------------------------------------------------------------------
// Constructor
pylith::problems::GreensFns::GreensFns(void) :
    _faultId(100),
    _faultImpulses(NULL),
    _integratorImpulses(NULL),
    _snes(NULL),
    _monitor(NULL),
    _residual(NULL),
    _solutionDot(NULL) {
    PyreComponent::setName(_GreensFns::pyreComponent);
} // constructor


// ---------------------------------------------------------------------------------------------------------------------
// Destructor
pylith::problems::GreensFns::~GreensFns(void) {
    deallocate();
} // destructor


// ---------------------------------------------------------------------------------------------------------------------
// Deallocate PETSc and local data structures.
void
pylith::problems::GreensFns::deallocate(void) {
    PYLITH_METHOD_BEGIN;

    Problem::deallocate();

    _faultImpulses = NULL; // Memory handle in Python. :TODO: Use shared pointer.
    _integratorImpulses = NULL; // Memory handle in Problem. :TODO: Use shared pointer.

    _monitor = NULL; // Memory handle in Python. :TODO: Use shared pointer.
    delete _residual;_residual = NULL;
    delete _solutionDot;_solutionDot = NULL;

    PetscErrorCode err = SNESDestroy(&_snes);PYLITH_CHECK_ERROR(err);

    PYLITH_METHOD_END;
} // deallocate


// ---------------------------------------------------------------------------------------------------------------------
// Set fault Id for problem.
void
pylith::problems::GreensFns::setFaultId(const int value) {
    PYLITH_COMPONENT_DEBUG("faultId(value="<<value<<")");

    _faultId = value;
} // setFaultId


// ---------------------------------------------------------------------------------------------------------------------
// Get fault Id for problem.
int
pylith::problems::GreensFns::getFaultId(void) const {
    return _faultId;
} // getFaultId


// ---------------------------------------------------------------------------------------------------------------------
// Set progress monitor.
void
pylith::problems::GreensFns::setProgressMonitor(pylith::problems::ProgressMonitorStep* monitor) {
    _monitor = monitor; // :KLUDGE: :TODO: Use shared pointer.
} // setProgressMonitor


// ---------------------------------------------------------------------------------------------------------------------
// Get Petsc DM associated with problem.
PetscDM
pylith::problems::GreensFns::getPetscDM(void) {
    PYLITH_METHOD_BEGIN;

    PetscDM dm = NULL;
    PetscErrorCode err = SNESGetDM(_snes, &dm);PYLITH_CHECK_ERROR(err);

    PYLITH_METHOD_RETURN(dm);
} // getPetscDM


// ---------------------------------------------------------------------------------------------------------------------
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
        if (_faultId == _interfaces[i]->getInterfaceId()) {
            faultImpulses = dynamic_cast<pylith::faults::FaultCohesiveImpulses*>(_interfaces[i]);
            if (!faultImpulses) {
                std::ostringstream msg;
                msg << "Found fault with id (" << _faultId << ") in interfaces for "
                    << "imposing impulses, but type is not FaultCohesiveImpulses.";
                throw std::runtime_error(msg.str());
            } // if
            break;
        } // if
    } // for
    if (!faultImpulses) {
        std::ostringstream msg;
        msg << "Could not find fault with id (" << _faultId << ") in interfaces for imposing impulses.";
        throw std::runtime_error(msg.str());
    } // if

    PYLITH_METHOD_END;
} // verifyConfiguration


// ---------------------------------------------------------------------------------------------------------------------
// Initialize.
void
pylith::problems::GreensFns::initialize(void) {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("initialize()");

    Problem::initialize();

    assert(_solution);

    // Find fault on which to put the impulses.
    assert(!_faultImpulses);
    const size_t numInterfaces = _interfaces.size();
    for (size_t i = 0; i < numInterfaces; ++i) {
        if (_faultId == _interfaces[i]->getInterfaceId()) {
            _faultImpulses = dynamic_cast<pylith::faults::FaultCohesiveImpulses*>(_interfaces[i]);
        } // if
    } // for
    assert(_faultImpulses);

    PetscErrorCode err = SNESDestroy(&_snes);PYLITH_CHECK_ERROR(err);assert(!_snes);
    const pylith::topology::Mesh& mesh = _solution->getMesh();
    err = SNESCreate(mesh.getComm(), &_snes);PYLITH_CHECK_ERROR(err);assert(_snes);
    err = SNESSetApplicationContext(_snes, (void*)this);PYLITH_CHECK_ERROR(err);
    err = SNESSetDM(_snes, _solution->getDM());PYLITH_CHECK_ERROR(err);

    PetscVec solutionVector = _solution->getGlobalVector();
    _solution->scatterLocalToVector(solutionVector);
    err = SNESSetSolution(_snes, solutionVector);PYLITH_CHECK_ERROR(err);
    delete _solutionDot;_solutionDot = new pylith::topology::Field(*_solution);assert(_solutionDot);
    _solutionDot->setLabel("solution_dot");

    // Initialize residual.
    delete _residual;_residual = new pylith::topology::Field(*_solution);assert(_residual);
    _residual->setLabel("residual");

    switch (_formulation) {
    case pylith::problems::Physics::QUASISTATIC:
        PYLITH_COMPONENT_DEBUG("Setting PetscSNES callbacks SNESSetFunction() and SNESSetJacobian().");
        err = SNESSetFunction(_snes, NULL, computeResidual, (void*)this);PYLITH_CHECK_ERROR(err);
        err = SNESSetJacobian(_snes, NULL, NULL, computeJacobian, (void*)this);PYLITH_CHECK_ERROR(err);
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

    err = SNESSetFromOptions(_snes);PYLITH_CHECK_ERROR(err);
    err = SNESSetUp(_snes);PYLITH_CHECK_ERROR(err);

    // Get integrator for fault with impulses.
    assert(!_integratorImpulses);
    const size_t numIntegrators = _integrators.size();
    for (size_t i = 0; i < numIntegrators; ++i) {
        if ((pylith::topology::Mesh::getCellsLabelName() == _integrators[i]->getLabelName()) &&
            ( _faultId == _integrators[i]->getLabelValue()) ) {
            _integratorImpulses = _integrators[i];
        } // if
    } // for
    if (!_integratorImpulses) {
        std::ostringstream msg;
        msg << "Could not find integrator for fault fault with id (" << _faultId << ") in integrators for problem.";
        throw std::runtime_error(msg.str());
    } // if

    pythia::journal::debug_t debug(pylith::utils::PyreComponent::getName());
    if (debug.state()) {
        PetscDS dsSoln = NULL;
        err = DMGetDS(_solution->getDM(), &dsSoln);PYLITH_CHECK_ERROR(err);
        debug << pythia::journal::at(__HERE__)
              << "Solution Discretization" << pythia::journal::endl;
        PetscDSView(dsSoln, PETSC_VIEWER_STDOUT_SELF);
    } // if

    if (_monitor) {
        _monitor->open();
    } // if

    PYLITH_METHOD_END;
} // initialize


// ---------------------------------------------------------------------------------------------------------------------
// Solve Green's functions problem.
void
pylith::problems::GreensFns::solve(void) {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("solve()");

    PetscErrorCode err;

    // Compute Jacobian
    PetscMat jacobian = NULL;
    PetscMat jacobianPrecond = NULL;
    createJacobian(jacobian, jacobianPrecond);
    err = SNESComputeJacobian(_snes, _solution->getGlobalVector(), jacobian, jacobianPrecond);PYLITH_CHECK_ERROR(err);

    PetscKSP ksp = NULL;
    err = SNESGetKSP(_snes, &ksp);PYLITH_CHECK_ERROR(err);assert(ksp);

    // Get number of impulses
    const size_t numImpulses = _faultImpulses->getNumImpulses();

    // Loop over impulses
    const PylithReal tolerance = 1.0e-4;
    for (size_t i = 0; i < numImpulses; ++i) {
        // Update impulse on fault
        const PetscReal impulseReal = i + tolerance;
        _integratorImpulses->updateState(impulseReal);

        err = KSPSolve(ksp, _residual->getGlobalVector(), _solution->getGlobalVector());PYLITH_CHECK_ERROR(err);
        _solution->scatterVectorToLocal(_solution->getGlobalVector());
        _solution->scatterLocalToOutput();
        poststep(impulseReal);
    } // for

    PYLITH_METHOD_END;
} // solve


// ---------------------------------------------------------------------------------------------------------------------
// Perform operations after advancing solution of one impulse.
void
pylith::problems::GreensFns::poststep(const double impulseReal) {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("poststep()");

    // Get current solution
    const PetscReal t = 0;
    const PetscReal dt = 1.0;

    // Update integrators.
    const size_t numIntegrators = _integrators.size();
    for (size_t i = 0; i < numIntegrators; ++i) {
        _integrators[i]->poststep(t, impulseReal, dt, *_solution);
    } // for

    // Update constraints.
    const size_t numConstraints = _constraints.size();
    for (size_t i = 0; i < numConstraints; ++i) {
        _constraints[i]->poststep(t, impulseReal, dt, *_solution);
    } // for

    // Notify problem observers of updated solution.
    assert(_observers);
    _observers->notifyObservers(t, (int)impulseReal, *_solution);

    // Get number of impulses for monitor
    const size_t numImpulses = _faultImpulses->getNumImpulses();
    if (_monitor) {
        assert(_normalizer);
        _monitor->update(impulseReal, 0.0, numImpulses*1.0);
    } // if

    PYLITH_METHOD_END;
} // poststep


// ----------------------------------------------------------------------
// Set solution values according to constraints (Dirichlet BC).
void
pylith::problems::GreensFns::setSolutionLocal(PetscVec solutionVec) {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("setSolutionLocal(solutionVec="<<solutionVec<<")");

    // Update PyLith view of the solution.
    assert(_solution);
    _solution->scatterVectorToLocal(solutionVec);

    const size_t numConstraints = _constraints.size();
    const PetscReal t = 0.0;
    for (size_t i = 0; i < numConstraints; ++i) {
        _constraints[i]->setSolution(_solution, t);
    } // for

    // _solution->view("SOLUTION AFTER SETTING VALUES");

    PYLITH_METHOD_END;
} // setSolutionLocal


// ----------------------------------------------------------------------
// Compute residual.
void
pylith::problems::GreensFns::computeResidual(PetscVec residualVec,
                                             PetscVec solutionVec) {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("computeResidual(solutionVec="<<solutionVec<<", residualVec="<<residualVec<<")");

    assert(residualVec);
    assert(solutionVec);
    assert(_solution);

    // Update PyLith view of the solution.
    setSolutionLocal(solutionVec);

    // Sum residual across integrators.
    const PetscReal t = 0.0;
    const PetscReal dt = 1.0;
    _residual->zeroLocal();
    const int numIntegrators = _integrators.size();
    assert(numIntegrators > 0); // must have at least 1 integrator
    for (int i = 0; i < numIntegrators; ++i) {
        _integrators[i]->computeLHSResidual(_residual, t, dt, *_solution, *_solutionDot);
    } // for

    // Assemble residual values across processes.
    PetscErrorCode err = VecSet(residualVec, 0.0);PYLITH_CHECK_ERROR(err);
    _residual->scatterLocalToVector(residualVec, ADD_VALUES);

    PYLITH_METHOD_END;
} // computeResidual


// ----------------------------------------------------------------------
// Compute Jacobian.
void
pylith::problems::GreensFns::computeJacobian(PetscMat jacobianMat,
                                             PetscMat precondMat,
                                             PetscVec solutionVec) {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("GreensFns::computeJacobian(solutionVec="<<solutionVec<<",jacobianMat="<<jacobianMat<<", precondMat="<<precondMat<<")");

    assert(jacobianMat);
    assert(precondMat);
    assert(solutionVec);

    // Zero Jacobian
    PetscErrorCode err = 0;
    PetscDS solnDS = NULL;
    PetscBool hasJacobian = PETSC_FALSE;
    err = DMGetDS(_solution->getDM(), &solnDS);PYLITH_CHECK_ERROR(err);
    err = PetscDSHasJacobian(solnDS, &hasJacobian);PYLITH_CHECK_ERROR(err);
    if (hasJacobian) { err = MatZeroEntries(jacobianMat);PYLITH_CHECK_ERROR(err); }
    err = MatZeroEntries(precondMat);PYLITH_CHECK_ERROR(err);

    // Update PyLith view of the solution.
    setSolutionLocal(solutionVec);

    // Sum Jacobian contributions across integrators.
    const PetscReal t = 0.0;
    const PetscReal dt = 1.0;
    const PetscReal s_tshift = 0.0;

    const size_t numIntegrators = _integrators.size();
    for (size_t i = 0; i < numIntegrators; ++i) {
        _integrators[i]->computeLHSJacobian(jacobianMat, precondMat, t, dt, s_tshift, *_solution, *_solutionDot);
    } // for

    // Solver handles assembly.

    PYLITH_METHOD_END;
} // computeJacobian


// ----------------------------------------------------------------------
// Create Jacobian. &jacobian, &jacobianPrecond
void
pylith::problems::GreensFns::createJacobian(PetscMat jacobianMat,
                                            PetscMat precondMat) {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("createJacobian(jacobianMat="<<jacobianMat<<", precondMat="<<precondMat<<")");

    assert(jacobianMat);
    assert(precondMat);

    PetscErrorCode err;
    PetscDS dsSoln = NULL;
    PetscBool hasJacobian, hasPreconditioner;
    PetscDM dmSNES = this->getPetscDM();
    err = DMCreateMatrix(dmSNES, &jacobianMat);PYLITH_CHECK_ERROR(err);
    err = DMGetDS(_solution->getDM(), &dsSoln);PYLITH_CHECK_ERROR(err);
    err = PetscDSHasJacobian(dsSoln, &hasJacobian);PYLITH_CHECK_ERROR(err);
    if (!hasJacobian) {
        throw std::runtime_error("PETSc DS indicates there is no Jacobian to form. Jacobian expected for Green's function problem.");
    } // if
    err = PetscDSHasJacobianPreconditioner(dsSoln, &hasPreconditioner);PYLITH_CHECK_ERROR(err);
    if (hasPreconditioner) {
        err = DMCreateMatrix(dmSNES, &precondMat);PYLITH_CHECK_ERROR(err);
        err = PetscObjectSetName((PetscObject) precondMat, "Preconditioning Matrix");PYLITH_CHECK_ERROR(err);
        err = PetscObjectSetOptionsPrefix((PetscObject) precondMat, "jacpre_");PYLITH_CHECK_ERROR(err);
    } // if

    PYLITH_METHOD_END;
} // createJacobian


// ---------------------------------------------------------------------------------------------------------------------
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


// ---------------------------------------------------------------------------------------------------------------------
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
