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

#include "Integrator.hh" // implementation of class methods

#include "pylith/topology/Mesh.hh" // USES Mesh
#include "pylith/topology/Field.hh" // USES Field
#include "pylith/problems/ObserversPhysics.hh" // USES ObserversPhysics
#include "pylith/problems/Physics.hh" // USES Physics

#include "pylith/utils/EventLogger.hh" // USES EventLogger
#include "pylith/utils/error.hh" // USES PYLITH_METHOD_*
#include "pylith/utils/journals.hh" // USES PYLITH_JOURNAL_*

#include "spatialdata/units/Nondimensional.hh" // USES Nondimensional

#include <cassert> // USES assert()
#include <typeinfo> // USES typeid()
#include <stdexcept> // USES std::runtime_error

// ---------------------------------------------------------------------------------------------------------------------
// Constructor
pylith::feassemble::Integrator::Integrator(pylith::problems::Physics* const physics) :
    PhysicsImplementation(physics),
    _labelName(""),
    _labelValue(1),
    _lhsJacobianTriggers(NEW_JACOBIAN_NEVER),
    _lhsJacobianLumpedTriggers(NEW_JACOBIAN_NEVER),
    _hasRHSResidual(false),
    _hasLHSResidual(false),
    _hasLHSJacobian(false),
    _hasLHSJacobianLumped(false),
    _needNewLHSJacobian(true),
    _needNewLHSJacobianLumped(true)
{}


// ---------------------------------------------------------------------------------------------------------------------
// Destructor
pylith::feassemble::Integrator::~Integrator(void) {
    deallocate();
} // destructor


// ---------------------------------------------------------------------------------------------------------------------
// Set name of label used to identify integration domain.
void
pylith::feassemble::Integrator::setLabelName(const char* name) {
    PYLITH_JOURNAL_DEBUG("setLabelName(name="<<name<<")");

    if (strlen(name) == 0) {
        throw std::runtime_error("Empty string given for name of label for integration domain.");
    } // if

    _labelName = name;
}


// ---------------------------------------------------------------------------------------------------------------------
// Get name of label used to identify integration domain.
const char*
pylith::feassemble::Integrator::getLabelName(void) const {
    return _labelName.c_str();
}


// ---------------------------------------------------------------------------------------------------------------------
// Set value of label used to identify integration domain.
void
pylith::feassemble::Integrator::setLabelValue(const int value) {
    PYLITH_JOURNAL_DEBUG("setLabelValue(value="<<value<<")");
    _labelValue = value;
}


// ---------------------------------------------------------------------------------------------------------------------
// Get value of label used to identify integration domain.
int
pylith::feassemble::Integrator::getLabelValue(void) const {
    return _labelValue;
}


// ---------------------------------------------------------------------------------------------------------------------
// Check whether LHS Jacobian needs to be recomputed.
bool
pylith::feassemble::Integrator::needNewLHSJacobian(const bool dtChanged) {
    if (_lhsJacobianTriggers & NEW_JACOBIAN_ALWAYS) {
        _needNewLHSJacobian = true;
    } else if (dtChanged && (_lhsJacobianTriggers & NEW_JACOBIAN_TIME_STEP_CHANGE)) {
        _needNewLHSJacobian = true;
    } // if

    return _needNewLHSJacobian;
} // needNewLHSJacobian


// ---------------------------------------------------------------------------------------------------------------------
// Check whether LHS lumped Jacobian needs to be recomputed.
bool
pylith::feassemble::Integrator::needNewLHSJacobianLumped(const bool dtChanged) {
    if (_lhsJacobianLumpedTriggers & NEW_JACOBIAN_ALWAYS) {
        _needNewLHSJacobianLumped = true;
    } else if (dtChanged && (_lhsJacobianLumpedTriggers & NEW_JACOBIAN_TIME_STEP_CHANGE)) {
        _needNewLHSJacobianLumped = true;
    } // if

    return _needNewLHSJacobianLumped;
} // needNewLHSJacobianLumped


// ---------------------------------------------------------------------------------------------------------------------
// Set LHS Jacobian trigger.
void
pylith::feassemble::Integrator::setLHSJacobianTriggers(const int value) {
    _lhsJacobianTriggers |= value;
} // setLHSJacobianTriggers


// ---------------------------------------------------------------------------------------------------------------------
// Set LHS lumped Jacobian trigger.
void
pylith::feassemble::Integrator::setLHSJacobianLumpedTriggers(const int value) {
    _lhsJacobianLumpedTriggers |= value;
} // setLHSJacobianLumpedTriggers


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


#include <iostream>

// ---------------------------------------------------------------------------------------------------------------------
// Update auxiliary field values to current time.
void
pylith::feassemble::Integrator::updateState(const PylithReal t) {
    PYLITH_METHOD_BEGIN;
    PYLITH_JOURNAL_DEBUG("updateState(t="<<t<<") empty method");

    PYLITH_METHOD_END;
} // updateState


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

    PetscDM dmSoln = solution.getDM();assert(dmSoln);
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
