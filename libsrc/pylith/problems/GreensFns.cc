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
    _faultImpulsesId(100),
    _faultImpulses(NULL),
    _snes(NULL),
    _monitor(NULL),
    _residual(NULL)
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

    _monitor = NULL; // Memory handle in Python. :TODO: Use shared pointer.
    delete _residual;_residual = NULL;

    PetscErrorCode err = SNESDestroy(&_snes);PYLITH_CHECK_ERROR(err);

    PYLITH_METHOD_END;
} // deallocate


// ---------------------------------------------------------------------------------------------------------------------
// Set fault Id for problem.
void
pylith::problems::GreensFns::setFaultId(const int value) {
    PYLITH_COMPONENT_DEBUG("setFaultId(value="<<value<<")");

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
/** Get nonlinear solver for problem.
 *
 * @returns PETSc SNES for problem.
 */
PetscSNES
pylith::problems::GreensFns::getPetscSNES(void) {
    return _snes;
} // getPetscSNES


// ---------------------------------------------------------------------------------------------------------------------
// Verify configuration.
void
pylith::problems::GreensFns::verifyConfiguration(void) const {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("Problem::verifyConfiguration(void)");

    Problem::verifyConfiguration();

    // Find fault on which to put the impulses.
    assert(!_faultImpulses);
    const size_t numInterfaces = _interfaces.size();
    for (size_t i = 0; i < numInterfaces; ++i) {
        if (_faultImpulsesId = _interfaces[i]->getInterfaceId()) {
            _faultImpulses = _interfaces[i];
        } // if
    } // for
    if (!_faultImpulses) {
        std::ostringstream msg;
        msg << "Could not find fault with id (" << _faultImpulsesId << ") in interfaces for imposing impulses.";
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

    PetscErrorCode err = SNESDestroy(&_snes);PYLITH_CHECK_ERROR(err);assert(!_snes);
    const pylith::topology::Mesh& mesh = _solution->mesh();
    err = SNESCreate(mesh.comm(), &_snes);PYLITH_CHECK_ERROR(err);assert(_snes);
    err = SNESSetApplicationContext(_snes, (void*)this);PYLITH_CHECK_ERROR(err);
    err = SNESSetDM(_snes, _solution->dmMesh());PYLITH_CHECK_ERROR(err);

    PetscVec solutionVector = _solution->globalVector();
    _solution->scatterLocalToVector(solutionVector);
    err = SNESSetSolution(_snes, solutionVector);PYLITH_CHECK_ERROR(err);

    // Initialize residual.
    delete _residual;_residual = new pylith::topology::Field(*_solution);assert(_residual);
    _residual->setLabel("residual");

    // Set callbacks.
    PYLITH_COMPONENT_DEBUG("Setting PetscSNES callback for poststep().");
    err = SNESSetPostStep(_snes, poststep);PYLITH_CHECK_ERROR(err);

    switch (_formulation) {
    case pylith::problems::Physics::QUASISTATIC:
        PYLITH_COMPONENT_DEBUG("Setting PetscSNES callbacks computeIFunction() and computeIJacobian().");
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

    pythia::journal::debug_t debug(pylith::utils::PyreComponent::getName());
    if (debug.state()) {
        PetscDS prob = NULL;
        err = DMGetDS(_solution->dmMesh(), &prob);PYLITH_CHECK_ERROR(err);
        debug << pythia::journal::at(__HERE__)
              << "Solution Discretization" << pythia::journal::endl;
        PetscDSView(prob, PETSC_VIEWER_STDOUT_SELF);
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
    err = SNESFormJacobian(_snes);PYLITH_CHECK_ERROR(err);

    PetscKSP ksp = NULL;
    err = SNESGetKSP(_snes, &ksp);PYLITH_CHECK_ERROR(err);assert(ksp);

    // Get number of impulses
    const size_t numImpulses = _faultImpulses->getNumImpulses();
    for (size_t i=0; i < numImpulses; ++i) {
        // Update impulse on fault
        _faultImpulses->update(i); // FIX THIS based on FaultCohesiveImpulses implementation

        err = KSPSolve(ksp, _residual->globalVector(), _solution->globalVector());PYLITH_CHECK_ERROR(err);
    _solution->scatterVectorToLocal(solution->globalVector());
    _solution->scatterLocalToOutput();
        poststep(i); 
    } // for

    PYLITH_METHOD_END;
} // solve


// ---------------------------------------------------------------------------------------------------------------------
// Perform operations after advancing solution of one impulse.
void
pylith::problems::GreensFns::poststep(const int impulse) {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("poststep()");

    // Get current solution
    PetscErrorCode err;
    PylithInt tindex = impulse;

    // Update integrators.
    const size_t numIntegrators = _integrators.size();
    for (size_t i = 0; i < numIntegrators; ++i) {
        _integrators[i]->poststep(tindex, *_solution);
    } // for

    // Update constraints.
    const size_t numConstraints = _constraints.size();
    for (size_t i = 0; i < numConstraints; ++i) {
        _constraints[i]->poststep(tindex, *_solution);
    } // for

    // Notify problem observers of updated solution.
    assert(_observers);
    _observers->notifyObservers(tindex, *_solution);

    if (_monitor) {
        assert(_normalizer);
        _monitor->update(tindex); // fix this
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
    for (size_t i = 0; i < numConstraints; ++i) {
        _constraints[i]->setSolution(_solution);
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
    _residual->zeroLocal();
    const int numIntegrators = _integrators.size();
    assert(numIntegrators > 0); // must have at least 1 integrator
    for (int i = 0; i < numIntegrators; ++i) {
        _integrators[i]->computeResidual(_residual, *_solution);
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
    PYLITH_COMPONENT_DEBUG("computeJacobian(solutionVec="<<solutionVec<<", jacobianMat="<<jacobianMat<<", precondMat="<<precondMat<<")");

    assert(jacobianMat);
    assert(precondMat);
    assert(solutionVec);

    // Zero Jacobian
    PetscErrorCode err = 0;
    PetscDS solnDS = NULL;
    PetscBool hasJacobian = PETSC_FALSE;
    err = DMGetDS(_solution->dmMesh(), &solnDS);PYLITH_CHECK_ERROR(err);
    err = PetscDSHasJacobian(solnDS, &hasJacobian);PYLITH_CHECK_ERROR(err);
    if (hasJacobian) { err = MatZeroEntries(jacobianMat);PYLITH_CHECK_ERROR(err); }
    err = MatZeroEntries(precondMat);PYLITH_CHECK_ERROR(err);

    // Update PyLith view of the solution.
    setSolutionLocal(solutionVec);

    // Sum Jacobian contributions across integrators.
    const size_t numIntegrators = _integrators.size();
    for (size_t i = 0; i < numIntegrators; ++i) {
        _integrators[i]->computeJacobian(jacobianMat, precondMat, *_solution);
    } // for

    // Solver handles assembly.

    PYLITH_METHOD_END;
} // computeJacobian


// ---------------------------------------------------------------------------------------------------------------------
// Callback static method for computing residual for LHS, F(t,s,\dot{s}).
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


// ---------------------------------------------------------------------------------------------------------------------
// Callback static method for operations after advancing solution of one impulse.
PetscErrorCode
pylith::problems::GreensFns::poststep(PetscSNES snes) {
    PYLITH_METHOD_BEGIN;
    pythia::journal::debug_t debug(_GreensFns::pyreComponent);
    debug << pythia::journal::at(__HERE__)
          << "poststep(snes="<<snes<<")" << pythia::journal::endl;

    GreensFns* problem = NULL;
    PetscErrorCode err = SNESGetApplicationContext(snes, (void*)&problem);PYLITH_CHECK_ERROR(err);assert(problem);
    problem->poststep();

    PYLITH_METHOD_RETURN(0);
} // poststep

// End of file
