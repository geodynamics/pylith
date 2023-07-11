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

#include "FieldTester.hh" // implementation of class methods

#include "pylith/topology/Field.hh" // USES Field
#include "pylith/topology/FieldQuery.hh" // USES FieldQuery
#include "pylith/utils/error.hh" // USES PYLITH_CHECK_ERROR

#include "spatialdata/spatialdb/SpatialDB.hh" // USES SpatialDB

#include "petscdm.h" // USES PetscDM

#include "catch2/catch_test_macros.hpp"

// ---------------------------------------------------------------------------------------------------------------------
// Check to make sure field matches spatial database.
PylithReal
pylith::testing::FieldTester::checkFieldWithDB(const pylith::topology::Field& field,
                                               spatialdata::spatialdb::SpatialDB* fieldDB,
                                               const PylithReal lengthScale) {
    PYLITH_METHOD_BEGIN;
    PylithReal norm = 0.0;
    PylithReal t = 0.0;

    const PetscDM dmField = field.getDM();assert(dmField);
    pylith::topology::FieldQuery fieldQuery(field);
    fieldQuery.initializeWithDefaultQueries();
    fieldQuery.openDB(fieldDB, lengthScale);
    PetscErrorCode err = DMPlexComputeL2DiffLocal(dmField, t, fieldQuery._functions, (void**)fieldQuery._contextPtrs,
                                                  field.getLocalVector(), &norm);assert(!err);
    fieldQuery.closeDB(fieldDB);

    PYLITH_METHOD_RETURN(norm);
} // checkFieldWithDB


// ---------------------------------------------------------------------------------------------------------------------
// Test subfield info created by factory.
void
pylith::testing::FieldTester::checkSubfieldInfo(const pylith::topology::Field& field,
                                                const pylith::topology::Field::SubfieldInfo& infoE) {
    PYLITH_METHOD_BEGIN;

    const pylith::topology::Field::SubfieldInfo& info = field.getSubfieldInfo(infoE.description.label.c_str());

    CHECK(infoE.index == info.index);
    INFO("Checking subfield " << infoE.description.label);

    // Description
    const pylith::topology::Field::Description& descriptionE = infoE.description;
    const pylith::topology::Field::Description& description = info.description;
    CHECK(descriptionE.label == description.label);
    CHECK(descriptionE.alias == description.alias);
    CHECK(descriptionE.vectorFieldType == description.vectorFieldType);
    CHECK(descriptionE.numComponents == description.numComponents);
    for (size_t i = 0; i < descriptionE.numComponents; ++i) {
        CHECK(descriptionE.componentNames[i] == description.componentNames[i]);
    } // for
    CHECK(descriptionE.scale == description.scale);
    CHECK(descriptionE.validator == description.validator);
    CHECK(descriptionE.hasHistory == description.hasHistory);
    CHECK(descriptionE.historySize == description.historySize);

    // Discretization
    const pylith::topology::Field::Discretization& feE = infoE.fe;
    const pylith::topology::Field::Discretization& fe = info.fe;
    CHECK(feE.basisOrder == fe.basisOrder);
    CHECK(feE.quadOrder == fe.quadOrder);
    CHECK(feE.dimension == fe.dimension);
    CHECK(feE.isFaultOnly == fe.isFaultOnly);
    CHECK(feE.cellBasis == fe.cellBasis);
    CHECK(feE.feSpace == fe.feSpace);
    CHECK(feE.isBasisContinuous == fe.isBasisContinuous);

    PYLITH_METHOD_END;
} // checkSubfieldInfo


// End of file
