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

#include "InitialConditionPatch.hh" // implementation of class methods

#include "pylith/topology/FieldQuery.hh" // USES FieldQuery
#include "pylith/topology/Field.hh" // USES Field
#include "pylith/topology/Mesh.hh" // USES Mesh

#include "spatialdata/spatialdb/SpatialDB.hh" // USES SpatialDB
#include "spatialdata/units/Nondimensional.hh" // USES Nondimensional

#include "pylith/utils/error.hh" // USES PYLITH_CHECK_ERROR
#include "pylith/utils/journals.hh" // USES PYLITH_COMPONENT_*
#include <cassert> // USES assert()

// ---------------------------------------------------------------------------------------------------------------------
// Constructor
pylith::problems::InitialConditionPatch::InitialConditionPatch(void) :
    _patchId(0),
    _db(NULL) {
    PyreComponent::setName("initialconditionpatch");
} // constructor


// ---------------------------------------------------------------------------------------------------------------------
// Destructor
pylith::problems::InitialConditionPatch::~InitialConditionPatch(void) {
    deallocate();
} // destructor


// ---------------------------------------------------------------------------------------------------------------------
// Deallocate PETSc and local data structures.
void
pylith::problems::InitialConditionPatch::deallocate(void) {
    PYLITH_METHOD_BEGIN;

    _db = NULL; // :KLLUDGE: Should use shared pointer.

    PYLITH_METHOD_END;
} // deallocate


// ---------------------------------------------------------------------------------------------------------------------
// Set material id associated with patch.
void
pylith::problems::InitialConditionPatch::setMaterialId(const int value) {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("setMaterialId(value="<<value<<")");

    _patchId = value;

    PYLITH_METHOD_END;
} // setMaterialId


// ---------------------------------------------------------------------------------------------------------------------
// Get material id associated with patch.
int
pylith::problems::InitialConditionPatch::getMaterialId(void) const {
    return _patchId;
} // getMaterialId


// ---------------------------------------------------------------------------------------------------------------------
// Set spatial database holding initial conditions.
void
pylith::problems::InitialConditionPatch::setDB(spatialdata::spatialdb::SpatialDB* db) {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("setDB(db="<<db<<")");

    _db = db;

    PYLITH_METHOD_END;
} // setDB


// ---------------------------------------------------------------------------------------------------------------------
// Verify configuration is acceptable.
void
pylith::problems::InitialConditionPatch::verifyConfiguration(const pylith::topology::Field& solution) const {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("verifyConfiguration(solution="<<solution.getLabel()<<")");

    InitialCondition::verifyConfiguration(solution);

    const PetscDM dmSoln = solution.getDM();
    PetscBool hasLabel = PETSC_FALSE;
    const char* const labelName = pylith::topology::Mesh::getCellsLabelName();
    PetscErrorCode err = DMHasLabel(dmSoln, labelName, &hasLabel);PYLITH_CHECK_ERROR(err);
    if (!hasLabel) {
        std::ostringstream msg;
        msg << "Could not find label '" << labelName << "' for setting patch for initial condition '"
            << PyreComponent::getIdentifier() << "'.";
        throw std::runtime_error(msg.str());
    } // if

    PetscDMLabel dmLabel = NULL;
    err = DMGetLabel(solution.getDM(), labelName, &dmLabel);PYLITH_CHECK_ERROR(err);
    PetscBool hasValue = PETSC_FALSE;
    err = DMLabelHasValue(dmLabel, _patchId, &hasValue);PYLITH_CHECK_ERROR(err);
    if (!hasValue) {
        std::ostringstream msg;
        msg << "Label '" << labelName << "' missing value '" << _patchId << "' for initial condition '"
            << PyreComponent::getIdentifier() << "'.";
        throw std::runtime_error(msg.str());
    } // if

    // DMLabelView(dmLabel, PETSC_VIEWER_STDOUT_SELF); // :DEBUG:

    PYLITH_METHOD_END;
} // verifyConfiguration


// ---------------------------------------------------------------------------------------------------------------------
// Set solver type.
void
pylith::problems::InitialConditionPatch::setValues(pylith::topology::Field* solution,
                                                   const spatialdata::units::Nondimensional& normalizer) {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("setValues(solution="<<solution<<", normalizer)");

    assert(solution);

    pylith::topology::FieldQuery fieldQuery(*solution);

    const size_t numSubfields = _subfields.size();
    for (size_t i = 0; i < numSubfields; ++i) {
        const char** queryValues = NULL;
        const size_t numValues = 0;
        const pylith::topology::FieldQuery::convertfn_type convertFn = NULL;
        fieldQuery.setQuery(_subfields[i].c_str(), queryValues, numValues, convertFn, _db);
    } // for

    fieldQuery.openDB(_db, normalizer.getLengthScale());
    const char* const labelName = pylith::topology::Mesh::getCellsLabelName();
    fieldQuery.queryDBLabel(labelName, _patchId);
    fieldQuery.closeDB(_db);

    pythia::journal::debug_t debug(PyreComponent::getName());
    if (debug.state()) {
        PYLITH_COMPONENT_DEBUG("Displaying solution field");
        solution->view("Solution field with initial values");
    } // if

    PYLITH_METHOD_END;
} // setValues


// End of file
