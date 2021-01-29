// -*- C++ -*-
//
// ----------------------------------------------------------------------
//
// Brad T. Aagaard, U.S. Geological Survey
// Charles A. Williams, GNS Science
// Matthew G. Knepley, University of Chicago
//
// This code was developed as part of the Computational Infrastructure
// for Geodynamics (http://geodynamics.org).
//
// Copyright (c) 2010-2016 University of California, Davis
//
// See COPYING for license information.
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
#include "pylith/utils/journals.hh" // USES PYLITH_JOURNAL_*

#include <cassert> // USES assert()
#include <typeinfo> // USES typeid()
#include <stdexcept> // USES std::runtime_error

// ---------------------------------------------------------------------------------------------------------------------
// Default constructor.
pylith::feassemble::ConstraintUserFn::ConstraintUserFn(pylith::problems::Physics* const physics) :
    Constraint(physics),
    _fn(NULL) {
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
    void* context = NULL;
    const PylithInt labelId = 1;
    err = DMGetDS(solution.dmMesh(), &prob);PYLITH_CHECK_ERROR(err);
    const PetscInt i_field = solution.subfieldInfo(_subfieldName.c_str()).index;
    err = PetscDSAddBoundary(prob, DM_BC_ESSENTIAL, _constraintLabel.c_str(), _constraintLabel.c_str(), i_field,
                             _constrainedDOF.size(), &_constrainedDOF[0], (void (*)(void))_fn, NULL, 1, &labelId, context);
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
    PetscDM dmSoln = solution->dmMesh();

    // Get label for constraint.
    PetscDMLabel dmLabel = NULL;
    err = DMGetLabel(dmSoln, _constraintLabel.c_str(), &dmLabel);PYLITH_CHECK_ERROR(err);

    void* context = NULL;
    const int labelId = 1;
    const int fieldIndex = solution->subfieldInfo(_subfieldName.c_str()).index;
    const PylithInt numConstrained = _constrainedDOF.size();
    assert(solution->localVector());
    err = DMPlexLabelAddCells(dmSoln, dmLabel);PYLITH_CHECK_ERROR(err);
    err = DMPlexInsertBoundaryValuesEssential(dmSoln, t, fieldIndex, numConstrained, &_constrainedDOF[0], dmLabel, 1,
                                              &labelId, _fn, context, solution->localVector());PYLITH_CHECK_ERROR(err);
    err = DMPlexLabelClearCells(dmSoln, dmLabel);PYLITH_CHECK_ERROR(err);

    pythia::journal::debug_t debug(GenericComponent::getName());
    if (debug.state()) {
        PYLITH_JOURNAL_DEBUG("Displaying solution field");
        solution->view("solution field");
    } // if

    PYLITH_METHOD_END;
} // setSolution


// End of file
