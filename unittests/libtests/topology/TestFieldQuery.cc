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
// Copyright (c) 2010-2017 University of California, Davis
//
// See COPYING for license information.
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

#include "spatialdata/geocoords/CSCart.hh" // USES CSCart
#include "spatialdata/spatialdb/UserFunctionDB.hh" // USES UserFunctionDB
#include "spatialdata/units/Nondimensional.hh" // USES Nondimensional


const double pylith::topology::TestFieldQuery::FILL_VALUE = -999.0;

// ----------------------------------------------------------------------
// Setup testing data.
void
pylith::topology::TestFieldQuery::setUp(void) {
    PYLITH_METHOD_BEGIN;

    _data = new TestFieldQuery_Data; CPPUNIT_ASSERT(_data);
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

    delete _data; _data = NULL;
    delete _mesh; _mesh = NULL;
    delete _field; _field = NULL;
    delete _query; _query = NULL;

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
// Test queryFn().
void
pylith::topology::TestFieldQuery::testQueryFn(void) {
    PYLITH_METHOD_BEGIN;

    _initialize();
    CPPUNIT_ASSERT(_query);

    const char* name = "abc";

    // Test with non-null query functions and database.
    pylith::topology::FieldQuery::queryfn_type fnA = NULL;
    spatialdata::spatialdb::UserFunctionDB dbA;
    _query->queryFn(name, fnA, &dbA);
    CPPUNIT_ASSERT_EQUAL(fnA, _query->queryFn(name));
    CPPUNIT_ASSERT_EQUAL((const spatialdata::spatialdb::SpatialDB*)&dbA, _query->queryDB(name));

    // Test with non-null query function but NULL database.
    pylith::topology::FieldQuery::queryfn_type fnB = NULL;
    _query->queryFn(name, fnB);
    CPPUNIT_ASSERT_EQUAL(fnB, _query->queryFn(name));
    CPPUNIT_ASSERT(!_query->queryDB(name));

    // Test with NULL query function and database.
    _query->queryFn(name, NULL);
    CPPUNIT_ASSERT(!_query->queryFn(name));
    CPPUNIT_ASSERT(!_query->queryDB(name));

    PYLITH_METHOD_END;
} // testQueryFn

// ----------------------------------------------------------------------
// Test openDB(), closeDB()..
void
pylith::topology::TestFieldQuery::testOpenClose(void) {
    PYLITH_METHOD_BEGIN;

    _initialize();
    CPPUNIT_ASSERT(_query);

    _query->initializeWithDefaultQueryFns();

    // Test with non-NULL database.
    _query->openDB(_data->auxDB, _data->normalizer->lengthScale());
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

    _query->initializeWithDefaultQueryFns();

    CPPUNIT_ASSERT(_data);
    CPPUNIT_ASSERT(_data->normalizer);
    _query->openDB(_data->auxDB, _data->normalizer->lengthScale());
    _query->queryDB();
    _query->closeDB(_data->auxDB);

    //_field->view("FIELD"); // :DEBUG:

    // Compute difference with respect to direct queries to database.
    // Unfortunately, this also uses a FieldQuery object.
    PylithReal norm = 0.0;
    const PylithReal t = 0.0;
    pylith::topology::FieldQuery query(*_field);
    query.initializeWithDefaultQueryFns();
    query.openDB(_data->auxDB, _data->normalizer->lengthScale());
    PetscErrorCode err = DMPlexComputeL2DiffLocal(_field->dmMesh(), t, query.functions(), (void**)query.contextPtrs(), _field->localVector(), &norm); CPPUNIT_ASSERT(!err);
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

    _query->initializeWithDefaultQueryFns();

    _query->openDB(NULL, 1.0);
    _query->queryDB();
    _query->closeDB(NULL);

    //_field->view("FIELD"); // :DEBUG:

    // Expect auxfield to still contain FILL_VALUE values.
    PetscErrorCode err;
    const PylithReal tolerance = 1.0e-6;

    PylithReal min = 0;
    err = VecMin(_field->localVector(), NULL, &min);CPPUNIT_ASSERT(!err);
    CPPUNIT_ASSERT_DOUBLES_EQUAL_MESSAGE("Auxiliary field values don't match FILL_VALUE value.", FILL_VALUE, min, abs(FILL_VALUE*tolerance));

    PylithReal max = 0.0;
    err = VecMax(_field->localVector(), NULL, &max);CPPUNIT_ASSERT(!err);
    CPPUNIT_ASSERT_DOUBLES_EQUAL_MESSAGE("Auxiliary field values don't match FILL_VALUE value.", FILL_VALUE, max, abs(FILL_VALUE*tolerance));

    PYLITH_METHOD_END;
} // testQuery

// ----------------------------------------------------------------------
// Test validatorPositive().
void
pylith::topology::TestFieldQuery::testValidatorPositive(void) {
    PYLITH_METHOD_BEGIN;

    CPPUNIT_ASSERT(0 < strlen(pylith::topology::FieldQuery::validatorPositive(-1.0)));
    CPPUNIT_ASSERT(0 < strlen(pylith::topology::FieldQuery::validatorPositive(0.0)));
    CPPUNIT_ASSERT(0 == strlen(pylith::topology::FieldQuery::validatorPositive(1.0)));

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
    }   // for

    size = numCells * numCorners;
    int_array cells(size);
    for (PylithInt i = 0; i < size; ++i) {
        cells[i] = _data->cells[i];
    }   // for

    delete _mesh; _mesh = new Mesh; CPPUNIT_ASSERT(_mesh);
    pylith::meshio::MeshBuilder::buildMesh(_mesh, &coordinates, numVertices, spaceDim, cells, numCells, numCorners, cellDim);

    CPPUNIT_ASSERT(_data->cs);
    _mesh->coordsys(_data->cs);
    CPPUNIT_ASSERT(_data->normalizer);
    pylith::topology::MeshOps::nondimensionalize(_mesh, *_data->normalizer);

    // Setup field
    delete _field; _field = new pylith::topology::Field(*_mesh);CPPUNIT_ASSERT(_field);
    _field->label("auxiliary test field");
    for (int i = 0; i < _data->numAuxSubfields; ++i) {
        CPPUNIT_ASSERT(_data->auxDescriptions);
        CPPUNIT_ASSERT(_data->auxDiscretizations);
        _field->subfieldAdd(_data->auxDescriptions[i], _data->auxDiscretizations[i]);
    } // for
    _field->subfieldsSetup();
    _field->allocate();
    PetscErrorCode err;
    err = VecSet(_field->localVector(), FILL_VALUE); CPPUNIT_ASSERT(!err);

    delete _query; _query = new pylith::topology::FieldQuery(*_field);

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
    auxDB(new spatialdata::spatialdb::UserFunctionDB)
{   // constructor
}   // constructor


// ----------------------------------------------------------------------
// Destructor
pylith::topology::TestFieldQuery_Data::~TestFieldQuery_Data(void) {
    cells = NULL; // Assigned from const.
    coordinates = NULL; // Assigned from const.

    delete cs; cs = NULL;
    delete normalizer; normalizer = NULL;

    auxSubfields = NULL; // Assigned from const.
    auxDescriptions = NULL; // Assigned from const.
    auxDiscretizations = NULL; // Assigned from const.
    delete auxDB; auxDB = NULL;
}   // destructor


// End of file
