// =================================================================================================
// This code is part of PyLith, developed through the Computational Infrastructure
// for Geodynamics (https://github.com/geodynamics/pylith).
//
// Copyright (c) 2010-2025, University of California, Davis and the PyLith Development Team.
// All rights reserved.
//
// See https://mit-license.org/ and LICENSE.md and for license information.
// =================================================================================================

#include <portinfo>

#include "pylith/feassemble/Integrator.hh" // implementation of class methods

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

// ------------------------------------------------------------------------------------------------
namespace pylith {
    namespace feassemble {
        class _Integrator {
public:

            class Events {
public:

                static
                void init(void);

                static pylith::utils::EventLogger logger;
                static PylithInt initialize;
                static PylithInt poststep;
                static PylithInt setState;
                static PylithInt updateStateVars;
                static PylithInt computeDiagnosticField;
                static PylithInt computeDerivedField;
            };

        }; // _Integrator
    } // feassemble
} // pylith

pylith::utils::EventLogger pylith::feassemble::_Integrator::Events::logger;
PylithInt pylith::feassemble::_Integrator::Events::initialize;
PylithInt pylith::feassemble::_Integrator::Events::poststep;

// ------------------------------------------------------------------------------------------------
void
pylith::feassemble::_Integrator::Events::init(void) {
    logger.setClassName("Integrator");
    logger.initialize();
    initialize = logger.registerEvent("PL:Integrator:initialize");
    poststep = logger.registerEvent("PL:Integrator:poststep");
}


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
    _needNewLHSJacobianLumped(true) {
    _Integrator::Events::init();
}


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
    PYLITH_JOURNAL_DEBUG("initialize(solution="<<solution.getLabel()<<")");
    _Integrator::Events::logger.eventBegin(_Integrator::Events::initialize);

    const pylith::topology::Mesh& physicsDomainMesh = getPhysicsDomainMesh();
    delete _auxiliaryField;_auxiliaryField = _physics->createAuxiliaryField(solution, physicsDomainMesh);
    delete _diagnosticField;_diagnosticField = _physics->createDiagnosticField(solution, physicsDomainMesh);
    _computeDiagnosticField();
    _observers = _physics->getObservers(); // Memory managed by Physics
    if (_observers) {
        _observers->setPhysicsImplementation(this);
        _observers->setTimeScale(_physics->getNormalizer().getTimeScale());

        const pylith::problems::Observer::NotificationType notification = pylith::problems::ObserverPhysics::DIAGNOSTIC;
        _observers->notifyObservers(0.0, 0, solution, notification);
    } // if
    delete _diagnosticField;_diagnosticField = NULL;

    delete _derivedField;_derivedField = _physics->createDerivedField(solution, physicsDomainMesh);

    _Integrator::Events::logger.eventEnd(_Integrator::Events::initialize);
    PYLITH_METHOD_END;
} // initialize


// ---------------------------------------------------------------------------------------------------------------------
// Set auxiliary field values for current time.
void
pylith::feassemble::Integrator::setState(const PylithReal t) {
    PYLITH_METHOD_BEGIN;
    PYLITH_JOURNAL_DEBUG("setState(t="<<t<<") empty method");

    PYLITH_METHOD_END;
} // setState


// ---------------------------------------------------------------------------------------------------------------------
// Update auxiliary fields at end of time step.
void
pylith::feassemble::Integrator::poststep(const PylithReal t,
                                         const PylithInt tindex,
                                         const PylithReal dt,
                                         const pylith::topology::Field& solution,
                                         const pylith::problems::Observer::NotificationType notification) {
    PYLITH_METHOD_BEGIN;
    PYLITH_JOURNAL_DEBUG("poststep(t="<<t<<", dt="<<dt<<")");
    _Integrator::Events::logger.eventBegin(_Integrator::Events::poststep);

    _updateStateVars(t, dt, solution);
    _computeDerivedField(t, dt, solution);
    notifyObservers(t, tindex, solution, notification);

    _Integrator::Events::logger.eventEnd(_Integrator::Events::poststep);
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
        PetscDS ds = NULL;
        err = DMGetRegionNumDS(dmSoln, i, label, fields, &ds, NULL);PYLITH_CHECK_ERROR(err);
        if (constants.size() > 0) {
            err = PetscDSSetConstants(ds, constants.size(), const_cast<double*>(&constants[0]));PYLITH_CHECK_ERROR(err);
        } else {
            err = PetscDSSetConstants(ds, 0, NULL);PYLITH_CHECK_ERROR(err);
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
// Compute diagnostic field from auxiliary field.
void
pylith::feassemble::Integrator::_computeDiagnosticField(void) {
    PYLITH_METHOD_BEGIN;
    PYLITH_JOURNAL_DEBUG("_computeDiagnosticField() empty method");

    // Default is to do nothing.

    PYLITH_METHOD_END;
} // _computeDiagnosticField


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
