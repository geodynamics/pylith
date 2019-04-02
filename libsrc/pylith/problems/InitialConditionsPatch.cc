// -*- C++ -*-
//
// ======================================================================
//
// Brad T. Aagaard, U.S. Geological Survey
// Charles A. Williams, GNS Science
// Matthew G. Knepley, University of Chicago
//
// This code was developed as part of the Computational Infrastructure
// for Geodynamics (http://geodynamics.org).
//
// Copyright (c) 2010-2015 University of California, Davis
//
// See COPYING for license information.
//
// ======================================================================
//

#include <portinfo>

#include "InitialConditionsPatch.hh" // implementation of class methods

#include "pylith/topology/FieldQuery.hh" // USES FieldQuery

#include "spatialdata/spatialdb/SpatialDB.hh" // USES SpatialDB
#include "spatialdata/units/Nondimensional.hh" // USES Nondimensional

#include "pylith/utils/error.hh" // USES PYLITH_CHECK_ERROR
#include "pylith/utils/journals.hh" // USES PYLITH_COMPONENT_*
#include <cassert> // USES assert()

// ----------------------------------------------------------------------
// Constructor
pylith::problems::InitialConditionsPatch::InitialConditionsPatch(void) :
    _patchLabel(""),
    _db(NULL) {
    PyreComponent::setName("initialconditionspatch");
} // constructor


// ---------------------------------------------------------------------------------------------------------------------
// Destructor
pylith::problems::InitialConditionsPatch::~InitialConditionsPatch(void) {
    deallocate();
} // destructor


// ---------------------------------------------------------------------------------------------------------------------
// Deallocate PETSc and local data structures.
void
pylith::problems::InitialConditionsPatch::deallocate(void) {
    PYLITH_METHOD_BEGIN;

    _db = NULL; // :KLLUDGE: Should use shared pointer.

    PYLITH_METHOD_END;
} // deallocate


// ---------------------------------------------------------------------------------------------------------------------
// Set label for marker associated with patch.
void
pylith::problems::InitialConditionsPatch::setMarkerLabel(const char* value) {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("setMarkerLabel(value="<<value<<")");

    if (strlen(value) == 0) {
        throw std::runtime_error("Empty string given for initial conditions patch label.");
    } // if

    _patchLabel = value;

    PYLITH_METHOD_END;
} // setMarkerLabel


// ---------------------------------------------------------------------------------------------------------------------
// Get label marking constrained degrees of freedom.
const char*
pylith::problems::InitialConditionsPatch::getMarkerLabel(void) const {
    return _patchLabel.c_str();
} // getMarkerLabel


// ---------------------------------------------------------------------------------------------------------------------
// Set spatial database holding initial conditions.
void
pylith::problems::InitialConditionsPatch::setDB(spatialdata::spatialdb::SpatialDB* db) {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("setDB(db="<<db<<")");

    _db = db;

    PYLITH_METHOD_END;
} // setDB


// ---------------------------------------------------------------------------------------------------------------------
// Set solver type.
void
pylith::problems::InitialConditionsPatch::setValues(pylith::topology::Field* solution,
                                                    const spatialdata::units::Nondimensional& normalizer) {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("setValues(solution="<<solution<<", normalizer)");

    assert(solution);

    pylith::topology::FieldQuery fieldQuery(*solution);
    fieldQuery.initializeWithDefaultQueryFns();
    fieldQuery.setMarkerLabel(_patchLabel.c_str());
    fieldQuery.openDB(_db, normalizer.lengthScale());
    fieldQuery.queryDB();
    fieldQuery.closeDB(_db);

    PYLITH_METHOD_END;
} // setValues


// End of file
