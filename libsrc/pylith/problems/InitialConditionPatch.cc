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
// Copyright (c) 2010-2022 University of California, Davis
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
#include "pylith/topology/Stratum.hh" // USES Stratum

#include "spatialdata/spatialdb/SpatialDB.hh" // USES SpatialDB
#include "spatialdata/units/Nondimensional.hh" // USES Nondimensional

#include "pylith/utils/error.hh" // USES PYLITH_CHECK_ERROR
#include "pylith/utils/journals.hh" // USES PYLITH_COMPONENT_*
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
    PYLITH_COMPONENT_DEBUG("setLabelName(value="<<value<<")");

    if (strlen(value) == 0) {
        throw std::runtime_error("Empty string given for label of initial condition patch.");
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
    PYLITH_COMPONENT_DEBUG("setDB(db="<<db<<")");

    _db = db;

    PYLITH_METHOD_END;
} // setDB


// ------------------------------------------------------------------------------------------------
// Verify configuration is acceptable.
void
pylith::problems::InitialConditionPatch::verifyConfiguration(const pylith::topology::Field& solution) const {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("verifyConfiguration(solution="<<solution.getLabel()<<")");

    InitialCondition::verifyConfiguration(solution);

    const PetscDM dmSoln = solution.getDM();
    PetscBool hasLabel = PETSC_FALSE;
    PetscErrorCode err = DMHasLabel(dmSoln, _labelName.c_str(), &hasLabel);PYLITH_CHECK_ERROR(err);
    if (!hasLabel) {
        std::ostringstream msg;
        msg << "Could not find label '" << _labelName << "' for setting patch for initial condition '"
            << PyreComponent::getIdentifier() << "'.";
        throw std::runtime_error(msg.str());
    } // if

    PetscDMLabel dmLabel = NULL;
    err = DMGetLabel(solution.getDM(), _labelName.c_str(), &dmLabel);PYLITH_CHECK_ERROR(err);assert(dmLabel);
    PetscBool hasValue = PETSC_FALSE;
    err = DMLabelHasValue(dmLabel, _labelValue, &hasValue);PYLITH_CHECK_ERROR(err);
    if (!hasValue) {
        std::ostringstream msg;
        msg << "Label '" << _labelName << "' missing value '" << _labelValue << "' for initial condition '"
            << PyreComponent::getIdentifier() << "'.";
        throw std::runtime_error(msg.str());
    } // if

    PetscInt stratumStart = -1, stratumEnd = -1;
    err = DMLabelGetStratumBounds(dmLabel, _labelValue, &stratumStart, &stratumEnd);PYLITH_CHECK_ERROR(err);
    pylith::topology::Stratum cellsStratum(dmSoln, pylith::topology::Stratum::HEIGHT, 0);
    if (stratumStart >= cellsStratum.begin() && stratumEnd <= cellsStratum.end()) {
        std::ostringstream msg;
        msg << "Label '" << _labelName << "' with value '" << _labelValue << "' for initial condition '"
            << PyreComponent::getIdentifier() << "' contains only cells. Labels for initial conditions must "
            "contain cells and lower dimension points (for example, vertices, edges, and faces). These are "
            "not available for CUBIT meshes; the are available for physical groups created using VertexGroup "
            "in Gmsh Python scripts.";
        throw std::runtime_error(msg.str());
    } // if

    PYLITH_METHOD_END;
} // verifyConfiguration


// ------------------------------------------------------------------------------------------------
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
    fieldQuery.queryDBLabel(_labelName.c_str(), _labelValue);
    fieldQuery.closeDB(_db);

    pythia::journal::debug_t debug(PyreComponent::getName());
    if (debug.state()) {
        PYLITH_COMPONENT_DEBUG("Displaying solution field");
        solution->view("Solution field with initial values");
    } // if

    PYLITH_METHOD_END;
} // setValues


// End of file
