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

#include "pylith/feassemble/UpdateStateVars.hh" // implementation of object methods

#include "pylith/topology/Field.hh" // USES Field

#include "pylith/utils/error.hh" // USES PYLITH_METHOD_*

#include <cassert> // USES assert()

// ---------------------------------------------------------------------------------------------------------------------
// Default constructor.
pylith::feassemble::UpdateStateVars::UpdateStateVars(void) :
    _stateVarsIS(NULL),
    _stateVarsDM(NULL),
    _stateVarsVecLocal(NULL),
    _stateVarsVecGlobal(NULL),
    _auxiliaryFieldVecGlobal(NULL) {}


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

    PylithCallPetsc(ISDestroy(&_stateVarsIS));
    PylithCallPetsc(DMDestroy(&_stateVarsDM));
    PylithCallPetsc(VecDestroy(&_stateVarsVecLocal));
    PylithCallPetsc(VecDestroy(&_stateVarsVecGlobal));
    PylithCallPetsc(VecDestroy(&_auxiliaryFieldVecGlobal));

    PYLITH_METHOD_END;
} // deallocate


// ---------------------------------------------------------------------------------------------------------------------
// Get PETSc DM associated with state variables.
PetscDM
pylith::feassemble::UpdateStateVars::stateVarsDM(void) {
    return _stateVarsDM;
} // stateVarsDM


// ---------------------------------------------------------------------------------------------------------------------
// Get PETSc local vector associated with state variables.
PetscVec
pylith::feassemble::UpdateStateVars::stateVarsLocalVector(void) {
    return _stateVarsVecLocal;
} // stateVarsLocalVector


// ---------------------------------------------------------------------------------------------------------------------
// Initialize layout for updating state variables.
void
pylith::feassemble::UpdateStateVars::initialize(const pylith::topology::Field& auxiliaryField) {
    PYLITH_METHOD_BEGIN;

    PetscDM auxiliaryDM = auxiliaryField.getDM();

    const pylith::string_vector& subfieldNames = auxiliaryField.getSubfieldNames();
    const size_t numAuxiliarySubfields = subfieldNames.size();
    pylith::int_array stateSubfieldIndices(numAuxiliarySubfields);

    size_t numStateSubfields = 0;
    for (size_t iSubfield = 0; iSubfield < numAuxiliarySubfields; ++iSubfield) {
        const pylith::topology::Field::SubfieldInfo& info = auxiliaryField.getSubfieldInfo(subfieldNames[iSubfield].c_str());
        if (info.description.hasHistory) {
            stateSubfieldIndices[numStateSubfields++] = info.index;
        } // if
    } // for
    std::sort(&stateSubfieldIndices[0], &stateSubfieldIndices[numStateSubfields]);

    // Create subDM holding only the state vars, which we want to update.
    PylithCallPetsc(DMCreateSubDM(auxiliaryDM, numStateSubfields, &stateSubfieldIndices[0], &_stateVarsIS,
                                  &_stateVarsDM));
    PylithCallPetsc(DMCreateGlobalVector(_stateVarsDM, &_stateVarsVecGlobal));
    PylithCallPetsc(DMCreateLocalVector(_stateVarsDM, &_stateVarsVecLocal));

    PylithCallPetsc(DMCreateGlobalVector(auxiliaryDM, &_auxiliaryFieldVecGlobal));

    PYLITH_METHOD_END;
} // initialize


// ---------------------------------------------------------------------------------------------------------------------
// Setup values for updating state variables.
void
pylith::feassemble::UpdateStateVars::prepare(pylith::topology::Field* auxiliaryField) {
    PYLITH_METHOD_BEGIN;

    // :TODO: Verify that we need the global vectors and can't get by with just using VecISCopy() with the local vector.

    PylithCallPetsc(VecSet(_stateVarsVecLocal, 0.0));

    // Move auxiliaryDM data to global vector.
    assert(auxiliaryField);
    PetscDM auxiliaryDM = auxiliaryField->getDM();
    PylithCallPetsc(DMLocalToGlobalBegin(auxiliaryDM, auxiliaryField->getLocalVector(), INSERT_VALUES, _auxiliaryFieldVecGlobal));
    PylithCallPetsc(DMLocalToGlobalEnd(auxiliaryDM, auxiliaryField->getLocalVector(), INSERT_VALUES, _auxiliaryFieldVecGlobal));

    PYLITH_METHOD_END;
} // prepare


// ---------------------------------------------------------------------------------------------------------------------
// Setup values for updating state variables.
void
pylith::feassemble::UpdateStateVars::restore(pylith::topology::Field* auxiliaryField) {
    PYLITH_METHOD_BEGIN;

    assert(auxiliaryField);
    PetscDM auxiliaryDM = auxiliaryField->getDM();

    // Move statevarDM data to global vector.
    PylithCallPetsc(DMLocalToGlobalBegin(_stateVarsDM, _stateVarsVecLocal, INSERT_VALUES, _stateVarsVecGlobal));
    PylithCallPetsc(DMLocalToGlobalEnd(_stateVarsDM, _stateVarsVecLocal, INSERT_VALUES, _stateVarsVecGlobal));

    // Copy global data from stateVars to auxiliaryField
    PylithCallPetsc(VecISCopy(_auxiliaryFieldVecGlobal, _stateVarsIS, SCATTER_FORWARD, _stateVarsVecGlobal));

    // Move auxiliaryDM data to local vector
    PylithCallPetsc(DMGlobalToLocalBegin(auxiliaryDM, _auxiliaryFieldVecGlobal, INSERT_VALUES, auxiliaryField->getLocalVector()));
    PylithCallPetsc(DMGlobalToLocalEnd(auxiliaryDM, _auxiliaryFieldVecGlobal, INSERT_VALUES, auxiliaryField->getLocalVector()));

    PYLITH_METHOD_END;
} // restore


// End of file
