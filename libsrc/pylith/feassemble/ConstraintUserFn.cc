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
#include "pylith/problems/ObserversPhysics.hh" // USES ObserversPhysics
#include "pylith/problems/Physics.hh" // USES Physics

#include "pylith/utils/EventLogger.hh" // USES EventLogger
#include "pylith/utils/error.hh" // USES PYLITH_METHOD_*
#include "pylith/utils/journals.hh" // USES PYLITH_JOURNAL_*

#include <cassert> // USES assert()
#include <typeinfo> // USES typeid()
#include <stdexcept> // USES std::runtime_error

// ---------------------------------------------------------------------------------------------------------------------
// Default constructor.
pylith::feassemble::ConstraintUserFn::ConstraintUserFn(pylith::problems::Physics* const physics) :
    Constraint(physics),
    _fn(NULL),
    _fnDot(NULL) {
    GenericComponent::setName("constraintuserfn");
} // constructor


// ---------------------------------------------------------------------------------------------------------------------
// Destructor.
pylith::feassemble::ConstraintUserFn::~ConstraintUserFn(void) {
    deallocate();
} // destructor


// ---------------------------------------------------------------------------------------------------------------------
// Set constraint kernel.
void
pylith::feassemble::ConstraintUserFn::setUserFn(const PetscUserFieldFunc fn) {
    _fn = fn;
} // setUserFn


// ---------------------------------------------------------------------------------------------------------------------
// Set constraint kernel time derivative.
void
pylith::feassemble::ConstraintUserFn::setUserFnDot(const PetscUserFieldFunc fnDot) {
    _fnDot = fnDot;
} // setUserFnDot


// ---------------------------------------------------------------------------------------------------------------------
// Initialize constraint domain, auxiliary field, and derived field. Update observers.
void
pylith::feassemble::ConstraintUserFn::initialize(const pylith::topology::Field& solution) {
    PYLITH_METHOD_BEGIN;
    PYLITH_JOURNAL_DEBUG("intialize(solution="<<solution.getLabel()<<")");

    Constraint::initialize(solution);

    const bool infoOnly = true;
    _observers->notifyObservers(0.0, 0, solution, infoOnly);

    // :KLUDGE: Potentially we may have multiple PetscDS objects. This assumes that the first one (with a NULL label) is
    // the correct one.
    PetscErrorCode err = 0;
    PetscDS prob = NULL;
    DMLabel label = NULL;
    void* context = NULL;
    const PylithInt labelId = 1;
    err = DMGetDS(solution.getDM(), &prob);PYLITH_CHECK_ERROR(err);
    const PetscInt i_field = solution.getSubfieldInfo(_subfieldName.c_str()).index;
    err = DMGetLabel(solution.getDM(), _constraintLabel.c_str(), &label);PYLITH_CHECK_ERROR(err);
    err = PetscDSAddBoundary(prob, DM_BC_ESSENTIAL, _constraintLabel.c_str(), label, 1, &labelId, i_field,
                             _constrainedDOF.size(), &_constrainedDOF[0], (void (*)(void))_fn, (void (*)(void))_fnDot, context, NULL);
    PYLITH_CHECK_ERROR(err);

    PYLITH_METHOD_END;
} // initialize


// ---------------------------------------------------------------------------------------------------------------------
// Set constrained values in solution field.
void
pylith::feassemble::ConstraintUserFn::setSolution(pylith::topology::Field* solution,
                                                  const double t) {
    PYLITH_METHOD_BEGIN;
    PYLITH_JOURNAL_DEBUG("setSolution(solution="<<solution->getLabel()<<", t="<<t<<")");

    assert(solution);

    PetscErrorCode err = 0;
    PetscDM dmSoln = solution->getDM();

    // Get label for constraint.
    PetscDMLabel dmLabel = NULL;
    err = DMGetLabel(dmSoln, _constraintLabel.c_str(), &dmLabel);PYLITH_CHECK_ERROR(err);

    void* context = NULL;
    const int labelId = 1;
    const int fieldIndex = solution->getSubfieldInfo(_subfieldName.c_str()).index;
    const PylithInt numConstrained = _constrainedDOF.size();
    assert(solution->getLocalVector());
    err = DMPlexLabelAddCells(dmSoln, dmLabel);PYLITH_CHECK_ERROR(err);
    err = DMPlexInsertBoundaryValuesEssential(dmSoln, t, fieldIndex, numConstrained, &_constrainedDOF[0], dmLabel, 1,
                                              &labelId, _fn, context, solution->getLocalVector());PYLITH_CHECK_ERROR(err);
    err = DMPlexLabelClearCells(dmSoln, dmLabel);PYLITH_CHECK_ERROR(err);

    pythia::journal::debug_t debug(GenericComponent::getName());
    if (debug.state()) {
        PYLITH_JOURNAL_DEBUG("Displaying solution field");
        solution->view("solution field");
    } // if

    PYLITH_METHOD_END;
} // setSolution


// ---------------------------------------------------------------------------------------------------------------------
// Set constrained values time derivative in solution field.
void
pylith::feassemble::ConstraintUserFn::setSolutionDot(pylith::topology::Field* solutionDot,
                                                     const double t) {
    PYLITH_METHOD_BEGIN;
    PYLITH_JOURNAL_DEBUG("setSolutionDot(solutionDot="<<solutionDot->getLabel()<<", t="<<t<<")");
    if (!_fnDot) { PYLITH_METHOD_END; }

    assert(solutionDot);

    PetscErrorCode err = 0;
    PetscDM dmSolnDot = solutionDot->getDM();

    // Get label for constraint.
    PetscDMLabel dmLabel = NULL;
    err = DMGetLabel(dmSolnDot, _constraintLabel.c_str(), &dmLabel);PYLITH_CHECK_ERROR(err);

    void* context = NULL;
    const int labelId = 1;
    const int fieldIndex = solutionDot->getSubfieldInfo(_subfieldName.c_str()).index;
    const PylithInt numConstrained = _constrainedDOF.size();
    assert(solutionDot->getLocalVector());
    err = DMPlexLabelAddCells(dmSolnDot, dmLabel);PYLITH_CHECK_ERROR(err);
    err = DMPlexInsertBoundaryValuesEssential(dmSolnDot, t, fieldIndex, numConstrained, &_constrainedDOF[0], dmLabel, 1,
                                              &labelId, _fnDot, context, solutionDot->getLocalVector());PYLITH_CHECK_ERROR(err);
    err = DMPlexLabelClearCells(dmSolnDot, dmLabel);PYLITH_CHECK_ERROR(err);

    pythia::journal::debug_t debug(GenericComponent::getName());
    if (debug.state()) {
        PYLITH_JOURNAL_DEBUG("Displaying solutionDot field");
        solutionDot->view("solutionDot field");
    } // if

    PYLITH_METHOD_END;
} // setSolution


// End of file
