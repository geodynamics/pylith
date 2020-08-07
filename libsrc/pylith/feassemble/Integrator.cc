// -*- C++ -*-
//
// ======================================================================
//
// Brad T. Aagaard, U.S. Geological Survey
// Charles A. Williams, GNS Science
// Matthew G. Knepley, Rice University
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

#include "Integrator.hh" // implementation of class methods

#include "pylith/topology/Mesh.hh" // USES Mesh
#include "pylith/topology/Field.hh" // USES Field
#include "pylith/problems/ObserversPhysics.hh" // USES ObserversPhysics
#include "pylith/problems/Physics.hh" // USES Physics

#include "pylith/utils/EventLogger.hh" // USES EventLogger
#include "pylith/utils/journals.hh" // USES PYLITH_JOURNAL_*

#include "spatialdata/units/Nondimensional.hh" // USES Nondimensional

#include <cassert> // USES assert()
#include <typeinfo> // USES typeid()
#include <stdexcept> // USES std::runtime_error

// ---------------------------------------------------------------------------------------------------------------------
// Constructor
pylith::feassemble::Integrator::Integrator(pylith::problems::Physics* const physics) :
    PhysicsImplementation(physics),
    _needNewRHSJacobian(true),
    _needNewLHSJacobian(true),
    _LHSJacobianTriggers(NEW_JACOBIAN_NEVER),
    _RHSJacobianTriggers(NEW_JACOBIAN_NEVER),
    _dtPrev(-1.0e+30)
{}


// ---------------------------------------------------------------------------------------------------------------------
// Destructor
pylith::feassemble::Integrator::~Integrator(void) {
    deallocate();
} // destructor


// ---------------------------------------------------------------------------------------------------------------------
// Check whether RHS Jacobian needs to be recomputed.
bool
pylith::feassemble::Integrator::needNewRHSJacobian(void) const {
    return _needNewRHSJacobian;
} // needNewRHSJacobian


// ---------------------------------------------------------------------------------------------------------------------
// Check whether LHS Jacobian needs to be recomputed.
bool
pylith::feassemble::Integrator::needNewLHSJacobian(void) const {
    return _needNewLHSJacobian;
} // needNewLHSJacobian


// ---------------------------------------------------------------------------------------------------------------------
// Set LHS Jacobian trigger.
void
pylith::feassemble::Integrator::setLHSJacobianTriggers(const NewJacobianTriggers value) {
    if (value == NEW_JACOBIAN_NEVER) {
        _LHSJacobianTriggers = value;
    } else {
        _LHSJacobianTriggers |= value;
    } // if/else
} // setLHSJacobianTriggers


// ---------------------------------------------------------------------------------------------------------------------
// Set RHS Jacobian trigger.
void
pylith::feassemble::Integrator::setRHSJacobianTriggers(const NewJacobianTriggers value) {
    if (value == NEW_JACOBIAN_NEVER) {
        _RHSJacobianTriggers = value;
    } else {
        _RHSJacobianTriggers |= value;
    } // if/else
} // setRHSJacobianTriggers


// ---------------------------------------------------------------------------------------------------------------------
// Initialize integration domain, auxiliary field, and derived field. Update observers.
void
pylith::feassemble::Integrator::initialize(const pylith::topology::Field& solution) {
    PYLITH_METHOD_BEGIN;
    PYLITH_JOURNAL_DEBUG("intialize(solution="<<solution.getLabel()<<")");

    const pylith::topology::Mesh& physicsDomainMesh = getPhysicsDomainMesh();

    delete _auxiliaryField;_auxiliaryField = _physics->createAuxiliaryField(solution, physicsDomainMesh);
    delete _derivedField;_derivedField = _physics->createDerivedField(solution, physicsDomainMesh);
    _observers = _physics->getObservers();assert(_observers); // Memory managed by Physics
    _observers->setPhysicsImplementation(this);

    const bool infoOnly = true;
    _observers->notifyObservers(0.0, 0, solution, infoOnly);
    _observers->setTimeScale(_physics->getNormalizer().getTimeScale());

    PYLITH_METHOD_END;
} // initialize


// ---------------------------------------------------------------------------------------------------------------------
// Update auxiliary fields at beginning of time step.
void
pylith::feassemble::Integrator::prestep(const PylithReal t,
                                        const PylithReal dt) {
    PYLITH_METHOD_BEGIN;
    PYLITH_JOURNAL_DEBUG("prestep(t="<<t<<", dt="<<dt<<")");

    if (dt != _dtPrev) {
        if (_LHSJacobianTriggers & NEW_JACOBIAN_TIME_STEP_CHANGE) {
            _needNewLHSJacobian = true;
        } // if
        if (_RHSJacobianTriggers & NEW_JACOBIAN_TIME_STEP_CHANGE) {
            _needNewRHSJacobian = true;
        } // if
        _dtPrev = dt;
    } // if

    PYLITH_METHOD_END;
} // prestep


// ---------------------------------------------------------------------------------------------------------------------
// Update auxiliary fields at end of time step.
void
pylith::feassemble::Integrator::poststep(const PylithReal t,
                                         const PylithInt tindex,
                                         const PylithReal dt,
                                         const pylith::topology::Field& solution) {
    PYLITH_METHOD_BEGIN;
    PYLITH_JOURNAL_DEBUG("poststep(t="<<t<<", dt="<<dt<<")");

    _updateStateVars(t, dt, solution);

    _computeDerivedField(t, dt, solution);
    notifyObservers(t, tindex, solution);

    PYLITH_METHOD_END;
} // poststep


// ---------------------------------------------------------------------------------------------------------------------
// Set constants used in finite-element kernels (point-wise functions).
void
pylith::feassemble::Integrator::_setKernelConstants(const pylith::topology::Field& solution,
                                                    const PylithReal dt) const {
    PYLITH_METHOD_BEGIN;
    PYLITH_JOURNAL_DEBUG("_setKernelConstants(solution="<<solution.getLabel()<<", dt="<<dt<<")");

    assert(_physics);
    const pylith::real_array& constants = _physics->getKernelConstants(dt);

    PetscDM dmSoln = solution.dmMesh();assert(dmSoln);
    PetscInt numDS = 0;
    PetscErrorCode err = DMGetNumDS(dmSoln, &numDS);PYLITH_CHECK_ERROR(err);
    for (PetscInt i = 0; i < numDS; ++i) {
        PetscDMLabel* label = NULL;
        PetscIS* fields = NULL;
        PetscDS prob = NULL;
        err = DMGetRegionNumDS(dmSoln, i, label, fields, &prob);PYLITH_CHECK_ERROR(err);
        if (constants.size() > 0) {
            err = PetscDSSetConstants(prob, constants.size(), const_cast<double*>(&constants[0]));PYLITH_CHECK_ERROR(err);
        } else {
            err = PetscDSSetConstants(prob, 0, NULL);PYLITH_CHECK_ERROR(err);
        } // if/else
    } // for

    PYLITH_METHOD_END;
} // _setKernelConstants


// ---------------------------------------------------------------------------------------------------------------------
// Update state variables as needed.
void
pylith::feassemble::Integrator::_updateStateVars(const PylithReal t,
                                                 const PylithReal dt,
                                                 const pylith::topology::Field& solution) {
    PYLITH_METHOD_BEGIN;
    PYLITH_JOURNAL_DEBUG("_updateStateVars(t="<<t<<", dt="<<dt<<", solution="<<solution.getLabel()<<") empty method");

    // Default is to do nothing.

    PYLITH_METHOD_END;
} // _updateStateVars


// ---------------------------------------------------------------------------------------------------------------------
// Compute field derived from solution and auxiliary field.
void
pylith::feassemble::Integrator::_computeDerivedField(const PylithReal t,
                                                     const PylithReal dt,
                                                     const pylith::topology::Field& solution) {
    PYLITH_METHOD_BEGIN;
    PYLITH_JOURNAL_DEBUG("_computeDerivedField(t="<<t<<", dt="<<dt<<", solution="<<solution.getLabel()<<") empty method");

    // Default is to do nothing.

    PYLITH_METHOD_END;
} // _computeDerivedField


// End of file
