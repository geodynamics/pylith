// -*- C++ -*-
//
// ----------------------------------------------------------------------
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
// ----------------------------------------------------------------------
//

#include <portinfo>

#include "pylith/feassemble/ConstraintUserFn.hh" // implementation of object methods

#include "pylith/topology/Mesh.hh" // USES Mesh
#include "pylith/topology/Field.hh" // USES Field
#include "pylith/feassemble/IntegrationData.hh" // USES IntegrationData
#include "pylith/problems/ObserversPhysics.hh" // USES ObserversPhysics
#include "pylith/problems/Physics.hh" // USES Physics

#include "pylith/utils/EventLogger.hh" // USES EventLogger
#include "pylith/utils/error.hh" // USES PYLITH_METHOD_*
#include "pylith/utils/journals.hh" // USES PYLITH_JOURNAL_*

#include <cassert> // USES assert()
#include <typeinfo> // USES typeid()
#include <stdexcept> // USES std::runtime_error

namespace pylith {
    namespace feassemble {
        class _ConstraintUserFn {
public:

            static
            void setSolution(const pylith::topology::Field* field,
                             const PylithReal t,
                             PetscUserFieldFunc fn,
                             const pylith::feassemble::ConstraintUserFn& constraint);

        };
    }
}

// ------------------------------------------------------------------------------------------------
// Default constructor.
pylith::feassemble::ConstraintUserFn::ConstraintUserFn(pylith::problems::Physics* const physics) :
    Constraint(physics),
    _fn(NULL),
    _fnDot(NULL) {
    GenericComponent::setName("constraintuserfn");
} // constructor


// ------------------------------------------------------------------------------------------------
// Destructor.
pylith::feassemble::ConstraintUserFn::~ConstraintUserFn(void) {
    deallocate();
} // destructor


// ------------------------------------------------------------------------------------------------
// Set constraint kernel.
void
pylith::feassemble::ConstraintUserFn::setUserFn(const PetscUserFieldFunc fn) {
    _fn = fn;
} // setUserFn


// ------------------------------------------------------------------------------------------------
// Set constraint kernel time derivative.
void
pylith::feassemble::ConstraintUserFn::setUserFnDot(const PetscUserFieldFunc fnDot) {
    _fnDot = fnDot;
} // setUserFnDot


// ------------------------------------------------------------------------------------------------
// Initialize constraint domain, auxiliary field, and derived field. Update observers.
void
pylith::feassemble::ConstraintUserFn::initialize(const pylith::topology::Field& solution) {
    PYLITH_METHOD_BEGIN;
    PYLITH_JOURNAL_DEBUG(_labelName<<"="<<_labelValue<<" intialize(solution="<<solution.getLabel()<<")");

    Constraint::initialize(solution);

    const bool infoOnly = true;
    _observers->notifyObservers(0.0, 0, solution, infoOnly);

    // :KLUDGE: Potentially we may have multiple PetscDS objects. This assumes that the first one (with a NULL
    // label) is the correct one.
    PetscErrorCode err = 0;
    PetscDS prob = NULL;
    DMLabel label = NULL;
    void* context = NULL;
    err = DMGetDS(solution.getDM(), &prob);PYLITH_CHECK_ERROR(err);
    const PetscInt i_field = solution.getSubfieldInfo(_subfieldName.c_str()).index;
    err = DMGetLabel(solution.getDM(), _labelName.c_str(), &label);PYLITH_CHECK_ERROR(err);
    err = PetscDSAddBoundary(prob, DM_BC_ESSENTIAL, _labelName.c_str(), label, 1, &_labelValue, i_field,
                             _constrainedDOF.size(), &_constrainedDOF[0], (void (*)(void)) _fn, (void (*)(void)) _fnDot, context, NULL);
    PYLITH_CHECK_ERROR(err);

    PYLITH_METHOD_END;
} // initialize


// ------------------------------------------------------------------------------------------------
// Set constrained values in solution field.
void
pylith::feassemble::ConstraintUserFn::setSolution(pylith::feassemble::IntegrationData* integrationData) {
    assert(integrationData);
    PYLITH_METHOD_BEGIN;
    PYLITH_JOURNAL_DEBUG(_labelName<<"="<<_labelValue<<" setSolution(integrationData="<<integrationData->str()<<")");

    const pylith::topology::Field* solution = integrationData->getField(pylith::feassemble::IntegrationData::solution);
    assert(solution);
    const PylithReal t = integrationData->getScalar(pylith::feassemble::IntegrationData::time);

    _ConstraintUserFn::setSolution(solution, t, _fn, *this);

    if (_fnDot && integrationData->hasField(pylith::feassemble::IntegrationData::solution_dot)) {
        const pylith::topology::Field* solutionDot = integrationData->getField(pylith::feassemble::IntegrationData::solution_dot);
        assert(solutionDot);
        _ConstraintUserFn::setSolution(solutionDot, t, _fnDot, *this);
    } // if

    PYLITH_METHOD_END;
} // setSolution


// ------------------------------------------------------------------------------------------------
// Set constrained values in solution field.
void
pylith::feassemble::_ConstraintUserFn::setSolution(const pylith::topology::Field* field,
                                                   const PylithReal t,
                                                   PetscUserFieldFunc fn,
                                                   const pylith::feassemble::ConstraintUserFn& constraint) {
    PYLITH_METHOD_BEGIN;
    assert(field);

    PetscErrorCode err = 0;
    PetscDM dmField = field->getDM();

    // Get label for constraint.
    PetscDMLabel dmLabel = NULL;
    err = DMGetLabel(dmField, constraint._labelName.c_str(), &dmLabel);PYLITH_CHECK_ERROR(err);

    void* context = NULL;
    const int fieldIndex = field->getSubfieldInfo(constraint._subfieldName.c_str()).index;
    const PylithInt numConstrained = constraint._constrainedDOF.size();
    assert(field->getLocalVector());
    err = DMPlexLabelAddCells(dmField, dmLabel);PYLITH_CHECK_ERROR(err);
    err = DMPlexInsertBoundaryValuesEssential(dmField, t, fieldIndex, numConstrained, &constraint._constrainedDOF[0], dmLabel, 1,
                                              &constraint._labelValue, fn, context, field->getLocalVector());PYLITH_CHECK_ERROR(err);
    err = DMPlexLabelClearCells(dmField, dmLabel);PYLITH_CHECK_ERROR(err);

    PYLITH_METHOD_END;
} // setSolution


// End of file
