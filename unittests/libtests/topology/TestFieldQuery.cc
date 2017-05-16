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
#include "pylith/topology/MeshOps.hh" // USES MeshOps::createDMMesh()
#include "pylith/topology/Stratum.hh" // USES Stratum
#include "pylith/topology/VisitorMesh.hh" // USES VecVisitorMesh

#include "pylith/meshio/MeshBuilder.hh" // Uses MeshBuilder
#include "pylith/utils/array.hh" // USES scalar_array

#include "spatialdata/geocoords/CSCart.hh" // USES CSCart

// ----------------------------------------------------------------------
// Setup testing data.
void
pylith::topology::TestFieldQuery::setUp(void)
{ // setUp
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
pylith::topology::TestFieldQuery::tearDown(void)
{ // tearDown
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
pylith::topology::TestFieldQuery::testConstructor(void)
{ // testConstructor
    PYLITH_METHOD_BEGIN;

    Mesh mesh;
    Field field(mesh);
    FieldQuery query(field);
    CPPUNIT_ASSERT(!query._functions);
    CPPUNIT_ASSERT(!query._contexts);
    CPPUNIT_ASSERT(!query._contextPtrs);

    PYLITH_METHOD_END;
} // testConstructor

// ----------------------------------------------------------------------
// Test queryFn().
void
pylith::topology::TestFieldQuery::testQueryFn(void)
{ // testQueryFn
    PYLITH_METHOD_BEGIN;

#if 0
    Mesh mesh;
    Field field(mesh);
    FieldQuery query(field);
    query.queryFn(_data->subfieldAName, _data->queryFnA);
    query.queryFn(_data->subfieldBName, _data->queryFnB);
    CPPUNIT_ASSERT(_data->queryFnB == query.queryFn(_data->subfieldBName));
    CPPUNIT_ASSERT(_data->queryFnA == query.queryFn(_data->subfieldAName));
#else
    CPPUNIT_ASSERT_MESSAGE(":TODO: @brad Test not implemented.", false);
#endif

    PYLITH_METHOD_END;
} // testQueryFn

// ----------------------------------------------------------------------
// Test openDB(), closeDB()..
void
pylith::topology::TestFieldQuery::testOpenClose(void)
{ // testOpenClose
    PYLITH_METHOD_BEGIN;

    CPPUNIT_ASSERT_MESSAGE(":TODO: @brad Test not implemented.", false);

    PYLITH_METHOD_END;
} // testOpenClose

// ----------------------------------------------------------------------
// Test queryDB().
void
pylith::topology::TestFieldQuery::testQuery(void)
{ // testQuery
    PYLITH_METHOD_BEGIN;

    CPPUNIT_ASSERT_MESSAGE(":TODO: @brad Test not implemented.", false);

    PYLITH_METHOD_END;
} // testQuery

// ----------------------------------------------------------------------
// Test dbQueryGeneric.
void
pylith::topology::TestFieldQuery::testDBQueryGeneric(void)
{ // testDBQueryGeneric
    PYLITH_METHOD_BEGIN;

    CPPUNIT_ASSERT_MESSAGE(":TODO: @brad Test not implemented.", false);

    PYLITH_METHOD_END;
} // testDBQueryGeneric

// ----------------------------------------------------------------------
// Test validatorPositive().
void
pylith::topology::TestFieldQuery::testValidatorPositive(void)
{ // testValidatorPositive
    PYLITH_METHOD_BEGIN;

    CPPUNIT_ASSERT_MESSAGE(":TODO: @brad Test not implemented.", false);

    PYLITH_METHOD_END;
} // testValidatorPositive

// ----------------------------------------------------------------------
void
pylith::topology::TestFieldQuery::_initialize(void)
{ // _initialize
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
    pylith::meshio::MeshBuilder::buildMesh(_mesh, &coordinates, numVertices, spaceDim, cells, numCells, numCorners,
                                           cellDim);

    spatialdata::geocoords::CSCart cs;
    cs.setSpaceDim(spaceDim);
    cs.initialize();
    _mesh->coordsys(&cs);

    // Setup field
    delete _field; _field = new Field(*_mesh);
    _field->label("solution");
    _field->subfieldAdd(_data->descriptionA, _data->discretizationA);
    _field->subfieldAdd(_data->descriptionB, _data->discretizationB);
    _field->subfieldsSetup();

    // Allocate field.
    _field->allocate();

    // Populate with values.
    PetscDM dmMesh = _mesh->dmMesh(); CPPUNIT_ASSERT(dmMesh);
    Stratum depthStratum(dmMesh, Stratum::DEPTH, 0);
    const PylithInt vStart = depthStratum.begin();
    const PylithInt vEnd = depthStratum.end();

    VecVisitorMesh fieldVisitor(*_field);
    const PylithInt fiberDim = _data->descriptionA.numComponents + _data->descriptionB.numComponents;
    PetscScalar* fieldArray = fieldVisitor.localArray();
    for (PylithInt v = vStart, indexA = 0, indexB = 0; v < vEnd; ++v) {
        // Set values for field A
        const PylithInt offA = fieldVisitor.sectionOffset(v);
        CPPUNIT_ASSERT_EQUAL(fiberDim, fieldVisitor.sectionDof(v));
        for (size_t d = 0; d < _data->descriptionA.numComponents; ++d) {
            fieldArray[offA+d] = _data->subfieldAValues[indexA++];
        } // for
          // Set values for field B
        const PylithInt offB = offA + _data->descriptionA.numComponents;
        for (size_t d = 0; d < _data->descriptionB.numComponents; ++d) {
            fieldArray[offB+d] = _data->subfieldBValues[indexB++];
        } // for
    } // for

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

    subfieldAValues(NULL),
    subfieldBValues(NULL)
{   // constructor
}   // constructor


// ----------------------------------------------------------------------
// Destructor
pylith::topology::TestFieldQuery_Data::~TestFieldQuery_Data(void)
{   // destructor
}   // destructor


// End of file
