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

#include "pylith/feassemble/ConstraintSimple.hh" // implementation of object methods

#include "spatialdata/units/Nondimensional.hh" // USES Nondimensional

#include "pylith/topology/Mesh.hh" // USES Mesh
#include "pylith/topology/Field.hh" // USES Field
#include "pylith/feassemble/IntegrationData.hh" // USES IntegrationData
#include "pylith/problems/ObserversPhysics.hh" // USES ObserversPhysics
#include "pylith/problems/Physics.hh" // USES Physics

#include "pylith/utils/EventLogger.hh" // USES EventLogger
#include "pylith/utils/error.hh" // USES PYLITH_METHOD_*
#include "pylith/utils/journals.hh" // USES PYLITH_JOURNAL_*

#include <cassert> // USES assert()
#include <stdexcept> // USES std::runtime_error

// ---------------------------------------------------------------------------------------------------------------------
// Default constructor.
pylith::feassemble::ConstraintSimple::ConstraintSimple(pylith::problems::Physics* const physics) :
    Constraint(physics),
    _fn(NULL) {
    GenericComponent::setName("constraintSimple");
} // constructor


// ---------------------------------------------------------------------------------------------------------------------
// Destructor.
pylith::feassemble::ConstraintSimple::~ConstraintSimple(void) {
    deallocate();
} // destructor


// ---------------------------------------------------------------------------------------------------------------------
// Set constraint kernel.
void
pylith::feassemble::ConstraintSimple::setUserFn(const PetscUserFieldFunc fn) {
    _fn = fn;
} // setSimple


// ---------------------------------------------------------------------------------------------------------------------
// Initialize constraint domain. Update observers.
void
pylith::feassemble::ConstraintSimple::initialize(const pylith::topology::Field& solution) {
    PYLITH_METHOD_BEGIN;
    PYLITH_JOURNAL_DEBUG("intialize(solution="<<solution.getLabel()<<")");

    assert(_physics);
    _observers = NULL;

    PetscErrorCode err = 0;
    PetscDM dm = solution.getDM();
    PetscDMLabel label;
    PetscDS ds = NULL;
    void* context = NULL;
    PetscInt i_field = -1;
    PetscInt *closure = NULL;
    PetscIS pointIS;
    const PetscInt *points;
    PetscInt point, cStart, cEnd, clSize;
    err = DMGetLabel(dm, _labelName.c_str(), &label);PYLITH_CHECK_ERROR(err);
    err = DMLabelGetStratumIS(label, _labelValue, &pointIS);PYLITH_CHECK_ERROR(err);
    err = ISGetIndices(pointIS, &points);PYLITH_CHECK_ERROR(err);
    point = points[0];
    err = ISRestoreIndices(pointIS, &points);PYLITH_CHECK_ERROR(err);
    err = ISDestroy(&pointIS);PYLITH_CHECK_ERROR(err);
    err = DMPlexGetHeightStratum(dm, 0, &cStart, &cEnd);PYLITH_CHECK_ERROR(err);
    err = DMPlexGetTransitiveClosure(dm, point, PETSC_FALSE, &clSize, &closure);PYLITH_CHECK_ERROR(err);
    for (int cl = 0; cl < clSize*2; cl += 2) {
        PetscDS cds;
        const PetscInt q = closure[cl];
        PetscInt Nf;

        if ((q < cStart) || (q >= cEnd)) { continue;}
        err = DMGetCellDS(dm, q, &cds);PYLITH_CHECK_ERROR(err);
        err = PetscDSGetNumFields(cds, &Nf);PYLITH_CHECK_ERROR(err);
        for (int f = 0; f < Nf; ++f) {
            PetscObject disc;
            const char *name;

            err = PetscDSGetDiscretization(cds, f, &disc);PYLITH_CHECK_ERROR(err);
            err = PetscObjectGetName(disc, &name);PYLITH_CHECK_ERROR(err);
            if (_subfieldName == std::string(name)) {ds = cds;i_field = f;break;}
        }
    }
    if (!ds) {
        std::ostringstream msg;

        msg << "INTERNAL ERROR in ConstraintSimple::initialize()\nCould not find a DS with a field named ''" << _subfieldName << "' in solution";
        throw std::logic_error(msg.str());
    }
    err = DMPlexRestoreTransitiveClosure(solution.getDM(), point, PETSC_FALSE, &clSize, &closure);PYLITH_CHECK_ERROR(err);
    err = DMGetLabel(solution.getDM(), _labelName.c_str(), &label);PYLITH_CHECK_ERROR(err);
    err = PetscDSAddBoundary(ds, DM_BC_ESSENTIAL, _labelName.c_str(), label, 1, &_labelValue, i_field,
                             _constrainedDOF.size(), &_constrainedDOF[0], (void (*)(void)) _fn, NULL, context, NULL);
    PYLITH_CHECK_ERROR(err);
    err = DMViewFromOptions(dm, NULL, "-constraint_simple_dm_view");PYLITH_CHECK_ERROR(err);
    {
        PetscInt Nds;
        err = DMGetNumDS(dm, &Nds);PYLITH_CHECK_ERROR(err);
        for (int s = 0; s < Nds; ++s) {
            err = DMGetRegionNumDS(dm, s, NULL, NULL, &ds);PYLITH_CHECK_ERROR(err);
            err = PetscObjectViewFromOptions((PetscObject) ds, NULL, "-constraint_simple_ds_view");PYLITH_CHECK_ERROR(err);
        }
    }

    PYLITH_METHOD_END;
} // initialize


// ---------------------------------------------------------------------------------------------------------------------
// Set constrained values in solution field.
void
pylith::feassemble::ConstraintSimple::setSolution(pylith::feassemble::IntegrationData* integrationData) {
    assert(integrationData);
    PYLITH_METHOD_BEGIN;
    PYLITH_JOURNAL_DEBUG("setSolution(integrationData="<<integrationData->str()<<")");

    const pylith::topology::Field* solution = integrationData->getField(pylith::feassemble::IntegrationData::solution);
    assert(solution);
    const PylithReal t = integrationData->getScalar(pylith::feassemble::IntegrationData::time);

    PetscErrorCode err = 0;
    PetscDM dmSoln = solution->getDM();

    // Get label for constraint.
    PetscDMLabel dmLabel = NULL;
    err = DMGetLabel(dmSoln, _labelName.c_str(), &dmLabel);PYLITH_CHECK_ERROR(err);

    void* context = NULL;
    const int _labelValue = 1;
    const int fieldIndex = solution->getSubfieldInfo(_subfieldName.c_str()).index;
    const PylithInt numConstrained = _constrainedDOF.size();
    assert(solution->getLocalVector());
    err = DMPlexLabelAddCells(dmSoln, dmLabel);PYLITH_CHECK_ERROR(err);
    err = DMPlexInsertBoundaryValuesEssential(dmSoln, t, fieldIndex, numConstrained, &_constrainedDOF[0], dmLabel, 1,
                                              &_labelValue, _fn, context, solution->getLocalVector());PYLITH_CHECK_ERROR(err);
    err = DMPlexLabelClearCells(dmSoln, dmLabel);PYLITH_CHECK_ERROR(err);

    pythia::journal::debug_t debug(GenericComponent::getName());
    if (debug.state()) {
        PYLITH_JOURNAL_DEBUG("Displaying solution field");
        solution->view("solution field");
    } // if

    PYLITH_METHOD_END;
} // setSolution


// End of file
