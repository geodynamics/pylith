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

#include "pylith/feassemble/ConstraintSimple.hh" // implementation of object methods

#include "pylith/topology/Mesh.hh" // USES Mesh
#include "pylith/topology/Field.hh" // USES Field
#include "pylith/feassemble/IntegrationData.hh" // USES IntegrationData
#include "pylith/problems/ObserversPhysics.hh" // USES ObserversPhysics
#include "pylith/problems/Physics.hh" // USES Physics

#include "pylith/utils/EventLogger.hh" // USES EventLogger
#include "pylith/utils/error.hh" // USES PYLITH_METHOD_*
#include "pylith/utils/journals.hh" // USES PYLITH_JOURNAL_*

#include "pylith/scales/Scales.hh" // USES Scales

#include <cassert> // USES assert()
#include <stdexcept> // USES std::runtime_error

// ---------------------------------------------------------------------------------------------------------------------
// Default constructor.
pylith::feassemble::ConstraintSimple::ConstraintSimple(pylith::problems::Physics* const physics) :
    Constraint(physics),
    _fn(NULL) {
    GenericComponent::setName("constraintSimple ");
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
    PYLITH_JOURNAL_DEBUG("initialize (solution = "<<solution.getLabel()<<") ");

    assert(_physics);
    _observers = NULL;

    PetscDM dm = solution.getDM();
    PetscDMLabel label;
    PetscDS ds = NULL;
    void* context = NULL;
    PetscInt i_field = -1;
    PetscInt *closure = NULL;
    PetscIS pointIS;
    const PetscInt *points;
    PetscInt point, cStart, cEnd, clSize;
    PylithCallPetsc(DMGetLabel(dm, _labelName.c_str(), &label));
    PylithCallPetsc(DMLabelGetStratumIS(label, _labelValue, &pointIS));
    if (!pointIS) {
        PYLITH_METHOD_END;
    } // if
    PylithCallPetsc(ISGetIndices(pointIS, &points));assert(points);
    point = points[0];
    PylithCallPetsc(ISRestoreIndices(pointIS, &points));
    PylithCallPetsc(ISDestroy(&pointIS));
    PylithCallPetsc(DMPlexGetHeightStratum(dm, 0, &cStart, &cEnd));
    PylithCallPetsc(DMPlexGetTransitiveClosure(dm, point, PETSC_FALSE, &clSize, &closure));
    for (int cl = 0; cl < clSize*2; cl += 2) {
        PetscDS cds;
        const PetscInt q = closure[cl];
        PetscInt Nf;

        if ((q < cStart) || (q >= cEnd)) { continue;}
        PylithCallPetsc(DMGetCellDS(dm, q, &cds, NULL));
        PylithCallPetsc(PetscDSGetNumFields(cds, &Nf));
        for (int f = 0; f < Nf; ++f) {
            PetscObject disc;
            const char *name;

            PylithCallPetsc(PetscDSGetDiscretization(cds, f, &disc));
            PylithCallPetsc(PetscObjectGetName(disc, &name));
            if (_subfieldName == std::string(name)) {ds = cds;i_field = f;break;}
        } // for
    } // for
    PetscInt numConstrainedDOF = _constrainedDOF.size();
    PetscInt* constrainedDOF = &_constrainedDOF[0];
    if (!ds) {
        // :KLUDGE: It is possible for a process to have a DOF that we need to constrain, but the process
        // may not have any cells with that DOF. The underlying code doesn't actually care if the point is
        // in the DS, so just get any DS and use it for the constraint.
        PylithCallPetsc(DMGetDS(dm, &ds));
        PetscInt numDS = 0, numFields = 0;
        i_field = solution.getSubfieldInfo(_subfieldName.c_str()).index;
        numConstrainedDOF = 0;
        constrainedDOF = NULL;

        PylithCallPetsc(DMGetNumDS(dm, &numDS));
        for (PetscInt s = 0; s < numDS; ++s) {
            PylithCallPetsc(DMGetRegionNumDS(dm, s, NULL, NULL, &ds, NULL));
            PylithCallPetsc(PetscDSGetNumFields(ds, &numFields));
            if (i_field < numFields) { break;}
        } // for

    } // if
    PylithCallPetsc(DMPlexRestoreTransitiveClosure(solution.getDM(), point, PETSC_FALSE, &clSize, &closure));
    PylithCallPetsc(DMGetLabel(solution.getDM(), _labelName.c_str(), &label));
    PylithCallPetsc(PetscDSAddBoundary(ds, DM_BC_ESSENTIAL, _labelName.c_str(), label, 1, &_labelValue, i_field,
                                       numConstrainedDOF, constrainedDOF, (void (*)(void)) _fn, NULL, context, NULL));
    PylithCallPetsc(DMViewFromOptions(dm, NULL, "-constraint_simple_dm_view "));
    {
        PetscInt numDS;
        PylithCallPetsc(DMGetNumDS(dm, &numDS));
        for (int s = 0; s < numDS; ++s) {
            PylithCallPetsc(DMGetRegionNumDS(dm, s, NULL, NULL, &ds, NULL));
            PylithCallPetsc(PetscObjectViewFromOptions((PetscObject) ds, NULL, "-constraint_simple_ds_view "));
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
    PYLITH_JOURNAL_DEBUG("setSolution (integrationData = "<<integrationData->str()<<") ");

    const pylith::topology::Field* solution = integrationData->getField(pylith::feassemble::IntegrationData::solution);
    assert(solution);
    const PylithReal t = integrationData->getScalar(pylith::feassemble::IntegrationData::time);

    PetscDM dmSoln = solution->getDM();

    // Get label for constraint.
    PetscDMLabel dmLabel = NULL;
    PylithCallPetsc(DMGetLabel(dmSoln, _labelName.c_str(), &dmLabel));

    void* context = NULL;
    const int _labelValue = 1;
    const int fieldIndex = solution->getSubfieldInfo(_subfieldName.c_str()).index;
    const PylithInt numConstrained = _constrainedDOF.size();
    assert(solution->getLocalVector());
    PylithCallPetsc(DMPlexLabelAddCells(dmSoln, dmLabel));
    PylithCallPetsc(DMPlexInsertBoundaryValuesEssential(dmSoln, t, fieldIndex, numConstrained, &_constrainedDOF[0], dmLabel, 1,
                                                        &_labelValue, _fn, context, solution->getLocalVector()));
    PylithCallPetsc(DMPlexLabelClearCells(dmSoln, dmLabel));

    pythia::journal::debug_t debug(GenericComponent::getName());
    if (debug.state()) {
        PYLITH_JOURNAL_DEBUG("Displaying solution field ");
        solution->view("solution field ");
    } // if

    PYLITH_METHOD_END;
} // setSolution


// End of file
