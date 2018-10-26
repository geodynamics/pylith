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

#include "FieldTester.hh" // implementation of class methods

#include "pylith/topology/Field.hh" // USES Field
#include "pylith/topology/FieldQuery.hh" // USES FieldQuery
#include "pylith/utils/error.hh" // USES PYLITH_CHECK_ERROR

#include "spatialdata/spatialdb/SpatialDB.hh" // USES SpatialDB

#include "petscdm.h" // USES PetscDM

#include <cppunit/extensions/HelperMacros.h>

// ---------------------------------------------------------------------------------------------------------------------
// Check to make sure field matches spatial database.
PylithReal
pylith::topology::FieldTester::checkFieldWithDB(const pylith::topology::Field& field,
                                                spatialdata::spatialdb::SpatialDB* fieldDB,
                                                const PylithReal lengthScale) {
    PYLITH_METHOD_BEGIN;
    PylithReal norm = 0.0;
    PylithReal t = 0.0;

    const PetscDM dmField = field.dmMesh();assert(dmField);
    pylith::topology::FieldQuery fieldQuery(field);
    fieldQuery.initializeWithDefaultQueryFns();
    fieldQuery.openDB(fieldDB, lengthScale);
    PetscErrorCode err = DMPlexComputeL2DiffLocal(dmField, t, fieldQuery.functions(), (void**)fieldQuery.contextPtrs(),
                                                  field.localVector(), &norm);CPPUNIT_ASSERT(!err);
    fieldQuery.closeDB(fieldDB);

    PYLITH_METHOD_RETURN(norm);
} // checkFieldWithDB


// ---------------------------------------------------------------------------------------------------------------------
// Test subfield info created by factory.
void
pylith::topology::FieldTester::checkSubfieldInfo(const pylith::topology::Field& field,
                                                 const pylith::topology::Field::SubfieldInfo& infoE) {
    PYLITH_METHOD_BEGIN;

    const pylith::topology::Field::SubfieldInfo& info = field.subfieldInfo(infoE.description.label.c_str());

    CPPUNIT_ASSERT_EQUAL(infoE.index, info.index);

    const std::string& msg = "Error checking subfield " + infoE.description.label;

    // Description
    const pylith::topology::Field::Description& descriptionE = infoE.description;
    const pylith::topology::Field::Description& description = info.description;
    CPPUNIT_ASSERT_EQUAL_MESSAGE(msg.c_str(), descriptionE.label, description.label);
    CPPUNIT_ASSERT_EQUAL_MESSAGE(msg.c_str(), descriptionE.alias, description.alias);
    CPPUNIT_ASSERT_EQUAL_MESSAGE(msg.c_str(), descriptionE.vectorFieldType, description.vectorFieldType);
    CPPUNIT_ASSERT_EQUAL_MESSAGE(msg.c_str(), descriptionE.numComponents, description.numComponents);
    for (size_t i = 0; i < descriptionE.numComponents; ++i) {
        CPPUNIT_ASSERT_EQUAL_MESSAGE(msg.c_str(), descriptionE.componentNames[i], description.componentNames[i]);
    } // for
    CPPUNIT_ASSERT_EQUAL_MESSAGE(msg.c_str(), descriptionE.scale, description.scale);
    CPPUNIT_ASSERT_EQUAL_MESSAGE(msg.c_str(), descriptionE.validator, description.validator);
    CPPUNIT_ASSERT_EQUAL_MESSAGE(msg.c_str(), descriptionE.hasHistory, description.hasHistory);
    CPPUNIT_ASSERT_EQUAL_MESSAGE(msg.c_str(), descriptionE.historySize, description.historySize);

    // Discretization
    const pylith::topology::Field::Discretization& feE = infoE.fe;
    const pylith::topology::Field::Discretization& fe = info.fe;
    CPPUNIT_ASSERT_EQUAL_MESSAGE(msg.c_str(), feE.basisOrder, fe.basisOrder);
    CPPUNIT_ASSERT_EQUAL_MESSAGE(msg.c_str(), feE.quadOrder, fe.quadOrder);
    CPPUNIT_ASSERT_EQUAL_MESSAGE(msg.c_str(), feE.isBasisContinuous, fe.isBasisContinuous);
    CPPUNIT_ASSERT_EQUAL_MESSAGE(msg.c_str(), feE.feSpace, fe.feSpace);

    PYLITH_METHOD_END;
} // checkSubfieldInfo


// End of file
