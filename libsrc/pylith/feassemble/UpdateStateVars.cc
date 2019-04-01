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

#include "UpdateStateVars.hh" // implementation of object methods

#include "pylith/topology/Field.hh" // USES Field

#include "pylith/utils/error.hh" // USES PYLITH_METHOD_*

#include <cassert> // USES assert()

// ---------------------------------------------------------------------------------------------------------------------
// Default constructor.
pylith::feassemble::UpdateStateVars::UpdateStateVars(void) :
    superIS(NULL),
    superDM(NULL),
    stateVarIS(NULL),
    stateVarDM(NULL),
    stateVarsSolnVecLocal(NULL),
    stateVarsSolnVecGlobal(NULL),
    stateVarsVecGlobal(NULL),
    auxiliaryFieldVecGlobal(NULL),
    solutionVecGlobal(NULL) {}


// ---------------------------------------------------------------------------------------------------------------------
// Destructor.
pylith::feassemble::UpdateStateVars::~UpdateStateVars(void) {
    deallocate();
} // destructor


// ---------------------------------------------------------------------------------------------------------------------
// Deallocate PETSc and local data structures.
void
pylith::feassemble::UpdateStateVars::deallocate(void) {
    PYLITH_METHOD_BEGIN;

    PetscErrorCode err = 0;
    err = ISDestroy(&superIS[0]);PYLITH_CHECK_ERROR(err);
    err = ISDestroy(&superIS[1]);PYLITH_CHECK_ERROR(err);
    err = DMDestroy(&superDM);PYLITH_CHECK_ERROR(err);
    err = ISDestroy(&stateVarIS);PYLITH_CHECK_ERROR(err);
    err = DMDestroy(&stateVarDM);PYLITH_CHECK_ERROR(err);
    err = VecDestroy(&stateVarsSolnVecLocal);PYLITH_CHECK_ERROR(err);
    err = VecDestroy(&stateVarsSolnVecGlobal);PYLITH_CHECK_ERROR(err);
    err = VecDestroy(&auxiliaryFieldVecGlobal);PYLITH_CHECK_ERROR(err);
    err = VecDestroy(&solutionVecGlobal);PYLITH_CHECK_ERROR(err);

    PYLITH_METHOD_END;
} // deallocate


// ---------------------------------------------------------------------------------------------------------------------
// Initialize layout for updating state variables.
void
pylith::feassemble::UpdateStateVars::initialize(const pylith::topology::Field& solution,
                                                const pylith::topology::Field& auxiliaryField) {
    PYLITH_METHOD_BEGIN;

    PetscErrorCode err = 0;
    PetscDM auxiliaryDM = auxiliaryField.dmMesh();
    PetscDM solutionDM = solution.dmMesh();

    const pylith::string_vector& subfieldNames = auxiliaryField.subfieldNames();
    const size_t numAuxiliarySubfields = subfieldNames.size();
    pylith::int_array stateSubfieldIndices(numAuxiliarySubfields);

    size_t numStateSubfields = 0;
    for (size_t iSubfield = 0; iSubfield < numAuxiliarySubfields; ++iSubfield) {
        const pylith::topology::Field::SubfieldInfo& info = auxiliaryField.subfieldInfo(subfieldNames[iSubfield].c_str());
        if (info.description.hasHistory) {
            stateSubfieldIndices[numStateSubfields++] = info.index;
        } // if
    } // for

    // Create subDM holding only the state vars, which we want to update.
    err = DMCreateSubDM(auxiliaryDM, numStateSubfields, &stateSubfieldIndices[0], &stateVarIS,
                        &stateVarDM);PYLITH_CHECK_ERROR(err);
    err = DMCreateGlobalVector(stateVarDM, &stateVarsVecGlobal);PYLITH_CHECK_ERROR(err);

    // Create superDM of {state vars, solution}
    PetscDM dms[2];
    dms[0] = stateVarDM;
    dms[1] = solutionDM;
    err = DMCreateSuperDM(dms, 2, &superIS, &superDM);PYLITH_CHECK_ERROR(err);
    err = DMCreateGlobalVector(superDM, &stateVarsSolnVecGlobal);PYLITH_CHECK_ERROR(err);
    err = DMCreateLocalVector(superDM, &stateVarsSolnVecLocal);PYLITH_CHECK_ERROR(err);

    // Copy state vars from auxiliary vars
    err = DMCreateGlobalVector(auxiliaryDM, &auxiliaryFieldVecGlobal);

    // Copy current state vars and solution into superDM space
    // these 2 lines in initialization.
    err = DMCreateGlobalVector(solutionDM, &solutionVecGlobal);

    PYLITH_METHOD_END;
} // initialize


// ---------------------------------------------------------------------------------------------------------------------
// Setup values for updating state variables.
void
pylith::feassemble::UpdateStateVars::prepareValues(pylith::topology::Field* auxiliaryField,
                                                   const pylith::topology::Field& solution) {
    PYLITH_METHOD_BEGIN;

    assert(auxiliaryField);

    PetscErrorCode err = 0;
    PetscDM auxiliaryDM = auxiliaryField->dmMesh();
    PetscDM solutionDM = solution.dmMesh();

    const pylith::string_vector& subfieldNames = auxiliaryField->subfieldNames();
    const size_t numAuxiliarySubfields = subfieldNames.size();
    pylith::int_array stateSubfieldIndices(numAuxiliarySubfields);

    err = DMLocalToGlobalBegin(auxiliaryDM, auxiliaryField->localVector(), INSERT_VALUES, auxiliaryFieldVecGlobal);PYLITH_CHECK_ERROR(err);
    err = DMLocalToGlobalEnd(auxiliaryDM, auxiliaryField->localVector(), INSERT_VALUES, auxiliaryFieldVecGlobal);PYLITH_CHECK_ERROR(err);
    err = VecISCopy(auxiliaryFieldVecGlobal, stateVarIS, SCATTER_REVERSE, stateVarsVecGlobal);PYLITH_CHECK_ERROR(err);

    err = DMLocalToGlobalBegin(solutionDM, solution.localVector(), INSERT_VALUES, solutionVecGlobal);PYLITH_CHECK_ERROR(err);
    err = DMLocalToGlobalEnd(solutionDM, solution.localVector(), INSERT_VALUES, solutionVecGlobal);PYLITH_CHECK_ERROR(err);
    err = VecISCopy(stateVarsSolnVecGlobal, superIS[0], SCATTER_FORWARD, stateVarsVecGlobal);PYLITH_CHECK_ERROR(err);
    err = VecISCopy(stateVarsSolnVecGlobal, superIS[1], SCATTER_FORWARD, solutionVecGlobal);PYLITH_CHECK_ERROR(err);

    // Move superDM data to a local vector
    err = DMGlobalToLocalBegin(superDM, stateVarsSolnVecGlobal, INSERT_VALUES, stateVarsSolnVecLocal);PYLITH_CHECK_ERROR(err);
    err = DMGlobalToLocalEnd(superDM, stateVarsSolnVecGlobal, INSERT_VALUES, stateVarsSolnVecLocal);PYLITH_CHECK_ERROR(err);

    // Attach superDM data as auxiliary for update state vars kernel
    err = PetscObjectCompose((PetscObject) auxiliaryDM, "dmAux", (PetscObject) superDM);PYLITH_CHECK_ERROR(err);
    err = PetscObjectCompose((PetscObject) auxiliaryDM, "A",     (PetscObject) stateVarsSolnVecLocal);PYLITH_CHECK_ERROR(err);

    PYLITH_METHOD_END;
} // prepareValues


// End of file
