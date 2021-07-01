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

#include "TestFieldQuery.hh" // Implementation of class methods

#include "pylith/topology/FieldQuery.hh" // USES FieldQuery
#include "pylith/topology/Field.hh" // USES Field
#include "pylith/topology/Mesh.hh" // USES Mesh

#include "pylith/meshio/MeshBuilder.hh" // Uses MeshBuilder
#include "pylith/topology/MeshOps.hh" // USES MeshOps
#include "pylith/utils/error.hh" // USES PYLITH_METHOD_*

#include "spatialdata/geocoords/CSCart.hh" // USES CSCart
#include "spatialdata/spatialdb/UserFunctionDB.hh" // USES UserFunctionDB
#include "spatialdata/units/Nondimensional.hh" // USES Nondimensional

const double pylith::topology::TestFieldQuery::FILL_VALUE = -999.0;

namespace pylith {
    namespace topology {
        namespace _TestFieldQuery {
            std::string
            converter(PylithScalar valueSubfield[],
                      const PylithInt numComponents,
                      const pylith::scalar_array dbValues,
                      const pylith::int_array dbIndices) {
                return std::string("Hello");
            }


        } // _TestFieldQuery
    } // topology
} // pylith

// ----------------------------------------------------------------------
// Setup testing data.
void
pylith::topology::TestFieldQuery::setUp(void) {
    PYLITH_METHOD_BEGIN;

    _data = new TestFieldQuery_Data;CPPUNIT_ASSERT(_data);
    _mesh = NULL;
    _field = NULL;
    _query = NULL;

    PYLITH_METHOD_END;
} // setUp


// ----------------------------------------------------------------------
// Tear down testing data.
void
pylith::topology::TestFieldQuery::tearDown(void) {
    PYLITH_METHOD_BEGIN;

    delete _data;_data = NULL;
    delete _mesh;_mesh = NULL;
    delete _field;_field = NULL;
    delete _query;_query = NULL;

    PYLITH_METHOD_END;
} // tearDown


// ----------------------------------------------------------------------
// Test constructor.
void
pylith::topology::TestFieldQuery::testConstructor(void) {
    PYLITH_METHOD_BEGIN;

    pylith::topology::Mesh mesh;
    pylith::topology::Field field(mesh);
    pylith::topology::FieldQuery query(field);
    CPPUNIT_ASSERT(!query._functions);
    CPPUNIT_ASSERT(!query._contexts);
    CPPUNIT_ASSERT(!query._contextPtrs);

    PYLITH_METHOD_END;
} // testConstructor


// ----------------------------------------------------------------------
// Test setQuery().
void
pylith::topology::TestFieldQuery::testSetQuery(void) {
    PYLITH_METHOD_BEGIN;

    _initialize();
    CPPUNIT_ASSERT(_query);

    { // Test with spatial database values, convert function, and database.
        const char* subfieldName = "ab";
        const size_t numDBValues = 2;
        const char* dbValues[numDBValues] = { "one", "two" };
        spatialdata::spatialdb::UserFunctionDB dbUser;
        _query->setQuery(subfieldName, dbValues, numDBValues, &_TestFieldQuery::converter, &dbUser);

        const pylith::topology::FieldQuery::SubfieldQuery& info = _query->_subfieldQueries[subfieldName];
        CPPUNIT_ASSERT_EQUAL_MESSAGE("Subfield ab: Mismatch in number of database values.",
                                     numDBValues, info.queryValues.size());
        for (size_t i = 0; i < numDBValues; ++i) {
            CPPUNIT_ASSERT_EQUAL_MESSAGE("Subfield ab: Mismatch in names of database values.",
                                         std::string(dbValues[i]), info.queryValues[i]);
        } // for
        CPPUNIT_ASSERT_EQUAL_MESSAGE("Subfield ab: Mismatch in convert function.",
                                     &_TestFieldQuery::converter, info.converter);
        CPPUNIT_ASSERT_EQUAL_MESSAGE("Subfield ab: Mismatch in database.",
                                     &dbUser, dynamic_cast<spatialdata::spatialdb::UserFunctionDB*>(info.db));
    }

    { // Test with spatial database values and database.
        const char* subfieldName = "cd";
        const size_t numDBValues = 3;
        const char* dbValues[numDBValues] = { "one", "two", "three" };
        spatialdata::spatialdb::UserFunctionDB dbUser;
        _query->setQuery(subfieldName, dbValues, numDBValues, NULL, &dbUser);

        const pylith::topology::FieldQuery::SubfieldQuery& info = _query->_subfieldQueries[subfieldName];
        CPPUNIT_ASSERT_EQUAL_MESSAGE("Subfield cd: Mismatch in number of database values.",
                                     numDBValues, info.queryValues.size());
        for (size_t i = 0; i < numDBValues; ++i) {
            CPPUNIT_ASSERT_EQUAL_MESSAGE("Subfield cd: Mismatch in names of database values.",
                                         std::string(dbValues[i]), info.queryValues[i]);
        } // for
        CPPUNIT_ASSERT_EQUAL_MESSAGE("Subfield cd: Mismatch in convert function.",
                                     pylith::topology::FieldQuery::convertfn_type(NULL), info.converter);
        CPPUNIT_ASSERT_EQUAL_MESSAGE("Subfield cd: Mismatch in database.",
                                     &dbUser, dynamic_cast<spatialdata::spatialdb::UserFunctionDB*>(info.db));
    }

    { // Test with spatial database values and convert function.
        const char* subfieldName = "ef";
        const size_t numDBValues = 1;
        const char* dbValues[numDBValues] = { "two" };
        _query->setQuery(subfieldName, dbValues, numDBValues, &_TestFieldQuery::converter);

        const pylith::topology::FieldQuery::SubfieldQuery& info = _query->_subfieldQueries[subfieldName];
        CPPUNIT_ASSERT_EQUAL_MESSAGE("Subfield ef: Mismatch in number of database values.",
                                     numDBValues, info.queryValues.size());
        for (size_t i = 0; i < numDBValues; ++i) {
            CPPUNIT_ASSERT_EQUAL_MESSAGE("Subfield ef: Mismatch in names of database values.",
                                         std::string(dbValues[i]), info.queryValues[i]);
        } // for
        CPPUNIT_ASSERT_EQUAL_MESSAGE("Subfield ef: Mismatch in convert function.",
                                     &_TestFieldQuery::converter, info.converter);
        CPPUNIT_ASSERT_EQUAL_MESSAGE("Subfield ef: Mismatch in database.",
                                     (spatialdata::spatialdb::SpatialDB*)NULL, info.db);
    }

    { // Test with spatial database values.
        const char* subfieldName = "gh";
        const size_t numDBValues = 3;
        const char* dbValues[numDBValues] = { "three", "two", "one" };
        _query->setQuery(subfieldName, dbValues, numDBValues, &_TestFieldQuery::converter);

        const pylith::topology::FieldQuery::SubfieldQuery& info = _query->_subfieldQueries[subfieldName];
        CPPUNIT_ASSERT_EQUAL_MESSAGE("Subfield gh: Mismatch in number of database values.",
                                     numDBValues, info.queryValues.size());
        for (size_t i = 0; i < numDBValues; ++i) {
            CPPUNIT_ASSERT_EQUAL_MESSAGE("Subfield gh: Mismatch in names of database values.",
                                         std::string(dbValues[i]), info.queryValues[i]);
        } // for
        CPPUNIT_ASSERT_EQUAL_MESSAGE("Subfield gh: Mismatch in convert function.",
                                     &_TestFieldQuery::converter, info.converter);
        CPPUNIT_ASSERT_EQUAL_MESSAGE("Subfield gh: Mismatch in database.",
                                     (spatialdata::spatialdb::SpatialDB*)NULL, info.db);
    }

    { // Test with defaults.
        const char* subfieldName = "displacement";
        const size_t numDBValuesE = 2;
        const char* dbValuesE[numDBValuesE] = { "displacement_x", "displacement_y" };

        _query->setQuery(subfieldName);

        const pylith::topology::FieldQuery::SubfieldQuery& info = _query->_subfieldQueries[subfieldName];
        CPPUNIT_ASSERT_EQUAL_MESSAGE("Subfield displacement: Mismatch in number of database values.",
                                     size_t(2), info.queryValues.size());
        for (size_t i = 0; i < numDBValuesE; ++i) {
            CPPUNIT_ASSERT_EQUAL_MESSAGE("Subfield displacement: Mismatch in names of database values.",
                                         std::string(dbValuesE[i]), info.queryValues[i]);
        } // for
        CPPUNIT_ASSERT_EQUAL_MESSAGE("Subfield displacement: Mismatch in convert function.",
                                     pylith::topology::FieldQuery::convertfn_type(NULL), info.converter);
        CPPUNIT_ASSERT_EQUAL_MESSAGE("Subfield displacement: Mismatch in database.",
                                     (spatialdata::spatialdb::SpatialDB*)NULL, info.db);
    }

    PYLITH_METHOD_END;
} // testSetQuery


// ----------------------------------------------------------------------
// Test openDB(), closeDB()..
void
pylith::topology::TestFieldQuery::testOpenClose(void) {
    PYLITH_METHOD_BEGIN;

    _initialize();
    CPPUNIT_ASSERT(_query);

    _query->initializeWithDefaultQueries();

    // Test with non-NULL database.
    _query->openDB(_data->auxDB, _data->normalizer->getLengthScale());
    _query->closeDB(_data->auxDB);
    // Nothing to verify.

    // Test with NULL database (should be okay).
    _query->openDB(NULL, 1.0);
    _query->closeDB(NULL);
    // Nothing to verify.

    PYLITH_METHOD_END;
} // testOpenClose


// ----------------------------------------------------------------------
// Test queryDB().
void
pylith::topology::TestFieldQuery::testQuery(void) {
    PYLITH_METHOD_BEGIN;

    _initialize();
    CPPUNIT_ASSERT(_query);
    CPPUNIT_ASSERT(_field);

    _query->initializeWithDefaultQueries();

    CPPUNIT_ASSERT(_data);
    CPPUNIT_ASSERT(_data->normalizer);
    _query->openDB(_data->auxDB, _data->normalizer->getLengthScale());
    _query->queryDB();
    _query->closeDB(_data->auxDB);

    // _field->view("FIELD"); // :DEBUG:

    // Compute difference with respect to direct queries to database.
    // Unfortunately, this also uses a FieldQuery object.
    PylithReal norm = 0.0;
    const PylithReal t = 0.0;
    pylith::topology::FieldQuery query(*_field);
    query.initializeWithDefaultQueries();
    query.openDB(_data->auxDB, _data->normalizer->getLengthScale());
    PetscErrorCode err = DMPlexComputeL2DiffLocal(_field->getDM(), t, query._functions, (void**)query._contextPtrs,
                                                  _field->getLocalVector(), &norm);CPPUNIT_ASSERT(!err);
    query.closeDB(_data->auxDB);
    const PylithReal tolerance = 1.0e-6;
    CPPUNIT_ASSERT_DOUBLES_EQUAL(0.0, norm, tolerance);

    PYLITH_METHOD_END;
} // testQuery


// ----------------------------------------------------------------------
// Test queryDB() with NULL database.
void
pylith::topology::TestFieldQuery::testQueryNull(void) {
    PYLITH_METHOD_BEGIN;

    _initialize();
    CPPUNIT_ASSERT(_query);

    _query->initializeWithDefaultQueries();

    _query->openDB(NULL, 1.0);
    _query->queryDB();
    _query->closeDB(NULL);

    // _field->view("FIELD"); // :DEBUG:

    // Expect auxfield to still contain FILL_VALUE values.
    PetscErrorCode err;
    const PylithReal tolerance = 1.0e-6;

    PylithReal min = 0;
    err = VecMin(_field->getLocalVector(), NULL, &min);CPPUNIT_ASSERT(!err);
    CPPUNIT_ASSERT_DOUBLES_EQUAL_MESSAGE("Auxiliary field values don't match FILL_VALUE value.", FILL_VALUE, min, abs(FILL_VALUE*tolerance));

    PylithReal max = 0.0;
    err = VecMax(_field->getLocalVector(), NULL, &max);CPPUNIT_ASSERT(!err);
    CPPUNIT_ASSERT_DOUBLES_EQUAL_MESSAGE("Auxiliary field values don't match FILL_VALUE value.", FILL_VALUE, max, abs(FILL_VALUE*tolerance));

    PYLITH_METHOD_END;
} // testQuery


// ----------------------------------------------------------------------
// Test validatorPositive().
void
pylith::topology::TestFieldQuery::testValidatorPositive(void) {
    PYLITH_METHOD_BEGIN;

    CPPUNIT_ASSERT(NULL != pylith::topology::FieldQuery::validatorPositive(-1.0));
    CPPUNIT_ASSERT(NULL != pylith::topology::FieldQuery::validatorPositive(0.0));
    CPPUNIT_ASSERT(NULL == pylith::topology::FieldQuery::validatorPositive(1.0));

    PYLITH_METHOD_END;
} // testValidatorPositive


// ----------------------------------------------------------------------
// Test validatorNonnegative().
void
pylith::topology::TestFieldQuery::testValidatorNonnegative(void) {
    PYLITH_METHOD_BEGIN;

    CPPUNIT_ASSERT(NULL != pylith::topology::FieldQuery::validatorNonnegative(-1.0));
    CPPUNIT_ASSERT(NULL == pylith::topology::FieldQuery::validatorNonnegative(0.0));
    CPPUNIT_ASSERT(NULL == pylith::topology::FieldQuery::validatorNonnegative(1.0));

    PYLITH_METHOD_END;
} // testValidatorPositive


// ----------------------------------------------------------------------
void
pylith::topology::TestFieldQuery::_initialize(void) {
    PYLITH_METHOD_BEGIN;
    CPPUNIT_ASSERT(_data);

    const int cellDim = _data->cellDim;
    const int numCells = _data->numCells;
    const int numVertices = _data->numVertices;
    const int numCorners = _data->numCorners;
    const int spaceDim = _data->cellDim;

    PylithInt size = numVertices * spaceDim;
    scalar_array coordinates(size);
    for (PylithInt i = 0; i < size; ++i) {
        coordinates[i] = _data->coordinates[i];
    } // for

    size = numCells * numCorners;
    int_array cells(size);
    for (PylithInt i = 0; i < size; ++i) {
        cells[i] = _data->cells[i];
    } // for

    delete _mesh;_mesh = new Mesh;CPPUNIT_ASSERT(_mesh);
    pylith::meshio::MeshBuilder::buildMesh(_mesh, &coordinates, numVertices, spaceDim, cells, numCells, numCorners, cellDim);

    CPPUNIT_ASSERT(_data->cs);
    _mesh->setCoordSys(_data->cs);
    CPPUNIT_ASSERT(_data->normalizer);
    pylith::topology::MeshOps::nondimensionalize(_mesh, *_data->normalizer);

    // Setup field
    delete _field;_field = new pylith::topology::Field(*_mesh);CPPUNIT_ASSERT(_field);
    _field->setLabel("auxiliary test field");
    for (int i = 0; i < _data->numAuxSubfields; ++i) {
        CPPUNIT_ASSERT(_data->auxDescriptions);
        CPPUNIT_ASSERT(_data->auxDiscretizations);
        _field->subfieldAdd(_data->auxDescriptions[i], _data->auxDiscretizations[i]);
    } // for
    _field->subfieldsSetup();
    _field->createDiscretization();
    _field->allocate();
    PetscErrorCode err;
    err = VecSet(_field->getLocalVector(), FILL_VALUE);CPPUNIT_ASSERT(!err);

    delete _query;_query = new pylith::topology::FieldQuery(*_field);

    PYLITH_METHOD_END;
} // _initialize


// ----------------------------------------------------------------------
// Constructor
pylith::topology::TestFieldQuery_Data::TestFieldQuery_Data(void) :
    cellDim(0),
    numVertices(0),
    numCells(0),
    numCorners(0),
    cells(NULL),
    coordinates(NULL),
    cs(new spatialdata::geocoords::CSCart),
    normalizer(new spatialdata::units::Nondimensional),
    numAuxSubfields(0),
    auxSubfields(NULL),
    auxDescriptions(NULL),
    auxDiscretizations(NULL),
    auxDB(new spatialdata::spatialdb::UserFunctionDB) { // constructor
} // constructor


// ----------------------------------------------------------------------
// Destructor
pylith::topology::TestFieldQuery_Data::~TestFieldQuery_Data(void) {
    cells = NULL; // Assigned from const.
    coordinates = NULL; // Assigned from const.

    delete cs;cs = NULL;
    delete normalizer;normalizer = NULL;

    auxSubfields = NULL; // Assigned from const.
    auxDescriptions = NULL; // Assigned from const.
    auxDiscretizations = NULL; // Assigned from const.
    delete auxDB;auxDB = NULL;
} // destructor


// End of file
