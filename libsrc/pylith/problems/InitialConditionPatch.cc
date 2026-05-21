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

#include "pylith/problems/InitialConditionPatch.hh" // implementation of class methods

#include "pylith/topology/FieldQuery.hh" // USES FieldQuery
#include "pylith/topology/Field.hh" // USES Field
#include "pylith/topology/Mesh.hh" // USES Mesh
#include "pylith/topology/Stratum.hh" // USES Stratum

#include "spatialdata/spatialdb/SpatialDB.hh" // USES SpatialDB
#include "pylith/scales/Scales.hh" // USES Scales

#include "pylith/utils/error.hh" // USES PylithCallPetsc()
#include "pylith/utils/journals.hh" // USES PYLITH_COMPONENT_*
#include "pylith/utils/Exceptions.hh" // USES Exception

#include <cassert> // USES assert()

// ------------------------------------------------------------------------------------------------
// Constructor
pylith::problems::InitialConditionPatch::InitialConditionPatch(void) :
    _labelName(pylith::topology::Mesh::cells_label_name),
    _labelValue(1),
    _db(NULL) {
    PyreComponent::setName("initialconditionpatch");
} // constructor


// ------------------------------------------------------------------------------------------------
// Destructor
pylith::problems::InitialConditionPatch::~InitialConditionPatch(void) {
    deallocate();
} // destructor


// ------------------------------------------------------------------------------------------------
// Deallocate PETSc and local data structures.
void
pylith::problems::InitialConditionPatch::deallocate(void) {
    PYLITH_METHOD_BEGIN;

    _db = NULL; // :KLLUDGE: Should use shared pointer.

    PYLITH_METHOD_END;
} // deallocate


// ------------------------------------------------------------------------------------------------
// Set name of label marking initial condition patch.
void
pylith::problems::InitialConditionPatch::setLabelName(const char* value) {
    PYLITH_COMPONENT_DEBUG(pylith::journal::application_flow, "setLabelName(value="<<value<<")");

    if (strlen(value) == 0) {
        PYLITH_COMPONENT_ERROR(pylith::ValueError, pylith::journal::user_input,
                               "Empty string given for label of initial condition patch.");
    } // if

    _labelName = value;
} // setLabelName


// ------------------------------------------------------------------------------------------------
// Get name of label marking initial condition patch.
const char*
pylith::problems::InitialConditionPatch::getLabelName(void) const {
    return _labelName.c_str();
} // getLabelName


// ------------------------------------------------------------------------------------------------
// Set value of label marking initial condition patch.
void
pylith::problems::InitialConditionPatch::setLabelValue(const int value) {
    _labelValue = value;
} // setLabelValue


// ------------------------------------------------------------------------------------------------
// Get value of label marking initial condition patch.
int
pylith::problems::InitialConditionPatch::getLabelValue(void) const {
    return _labelValue;
} // getLabelValue


// ------------------------------------------------------------------------------------------------
// Set spatial database holding initial conditions.
void
pylith::problems::InitialConditionPatch::setDB(spatialdata::spatialdb::SpatialDB* db) {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG(pylith::journal::application_flow, "setDB(db="<<db<<")");

    _db = db;

    PYLITH_METHOD_END;
} // setDB


// ------------------------------------------------------------------------------------------------
// Verify configuration is acceptable.
void
pylith::problems::InitialConditionPatch::verifyConfiguration(const pylith::topology::Field& solution) const {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG(pylith::journal::application_flow, "verifyConfiguration(solution="<<solution.getLabel()<<")");

    InitialCondition::verifyConfiguration(solution);

    const PetscDM dmSoln = solution.getDM();
    PetscBool hasLabel = PETSC_FALSE;
    PylithCallPetsc(DMHasLabel(dmSoln, _labelName.c_str(), &hasLabel));
    if (!hasLabel) {
        PYLITH_COMPONENT_ERROR(pylith::ValueError, pylith::journal::user_input,
                               "Could not find label '" << _labelName << "' for setting patch for initial condition.");
    } // if

    PetscDMLabel dmLabel = NULL;
    PylithCallPetsc(DMGetLabel(solution.getDM(), _labelName.c_str(), &dmLabel));assert(dmLabel);
    PetscBool hasLabelValue = PETSC_FALSE;
    PylithCallPetsc(DMLabelHasValue(dmLabel, _labelValue, &hasLabelValue));
    int hasLabelValueIntLocal = int(hasLabelValue);
    int hasLabelValueInt = 0;
    PylithCallPetsc(MPI_Allreduce(&hasLabelValueIntLocal, &hasLabelValueInt, 1, MPI_INT, MPI_MAX,
                                  PetscObjectComm((PetscObject) dmSoln)));
    if (!hasLabelValueInt) {
        PYLITH_COMPONENT_ERROR(pylith::ValueError, pylith::journal::user_input,
                               "Label '" << _labelName << "' missing value '" << _labelValue << "' for initial condition.");
    } // if

    PetscInt stratumStart = -1, stratumEnd = -1;
    PylithCallPetsc(DMLabelGetStratumBounds(dmLabel, _labelValue, &stratumStart, &stratumEnd));
    pylith::topology::Stratum cellsStratum(dmSoln, pylith::topology::Stratum::HEIGHT, 0);
    if ((stratumStart >= cellsStratum.begin()) && (stratumEnd <= cellsStratum.end())) {
        PYLITH_COMPONENT_ERROR(pylith::ValueError, pylith::journal::user_input,
                               "Label '" << _labelName << "' with value '" << _labelValue << "' for initial condition contains only cells. "
                                         << "Labels for initial conditions must contain cells and lower dimension points "
                                         << "(for example, vertices, edges, and faces). These are not available for Cubit meshes; "
                                         << "they are available for physical groups created in Gmsh.");
    } // if

    PYLITH_METHOD_END;
} // verifyConfiguration


// ------------------------------------------------------------------------------------------------
// Set solver type.
void
pylith::problems::InitialConditionPatch::setValues(pylith::topology::Field* solution,
                                                   const pylith::scales::Scales& scales) {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG(pylith::journal::application_flow, "setValues(solution="<<solution<<", scales)");

    assert(solution);

    pylith::topology::FieldQuery fieldQuery(*solution);

    const size_t numSubfields = _subfields.size();
    for (size_t i = 0; i < numSubfields; ++i) {
        const char** queryValues = NULL;
        const size_t numValues = 0;
        const pylith::topology::FieldQuery::convertfn_type convertFn = NULL;
        fieldQuery.setQuery(_subfields[i].c_str(), queryValues, numValues, convertFn, _db);
    } // for

    fieldQuery.openDB(_db, scales.getLengthScale());
    fieldQuery.queryDBLabel(_labelName.c_str(), _labelValue);
    fieldQuery.closeDB(_db);

    pythia::journal::debug_t debug(pylith::journal::solution);
    if (debug.state()) {
        PYLITH_COMPONENT_DEBUG(pylith::journal::solution, "Displaying solution field");
        solution->view("Solution field with initial condition values");
    } // if

    PYLITH_METHOD_END;
} // setValues


// End of file
