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

#include "TestFieldsSubmesh.hh" // Implementation of class methods

#include "pylith/topology/Mesh.hh" // USES Mesh
#include "pylith/topology/MeshOps.hh" // USES createLowerDimMesh()
#include "pylith/topology/Field.hh" // USES Field
#include "pylith/topology/Fields.hh" // USES Fields
#include "pylith/topology/Stratum.hh" // USES Stratum
#include "pylith/topology/VisitorMesh.hh" // USES VecVisitorMesh

#include "pylith/meshio/MeshIOAscii.hh" // USES MeshIOAscii

// ----------------------------------------------------------------------
CPPUNIT_TEST_SUITE_REGISTRATION(pylith::topology::TestFieldsSubmesh);

// ----------------------------------------------------------------------
void
pylith::topology::TestFieldsSubmesh::setUp(void) {
    PYLITH_METHOD_BEGIN;

    _mesh = new Mesh();
    meshio::MeshIOAscii importer;
    importer.filename("data/tri3.mesh");
    importer.read(_mesh);

    _submesh = MeshOps::createLowerDimMesh(*_mesh, "bc");CPPUNIT_ASSERT(_submesh);

    PYLITH_METHOD_END;
} // setUp


// ----------------------------------------------------------------------
void
pylith::topology::TestFieldsSubmesh::tearDown(void) {
    PYLITH_METHOD_BEGIN;

    delete _mesh;_mesh = NULL;
    delete _submesh;_submesh = NULL;

    PYLITH_METHOD_END;
} // tearDown


// ----------------------------------------------------------------------
// Test constructor.
void
pylith::topology::TestFieldsSubmesh::testConstructor(void) {
    PYLITH_METHOD_BEGIN;

    CPPUNIT_ASSERT(_submesh);
    Fields fields(*_submesh);

    PYLITH_METHOD_END;
} // testConstructor


// ----------------------------------------------------------------------
// Test add().
void
pylith::topology::TestFieldsSubmesh::testAdd(void) {
    PYLITH_METHOD_BEGIN;

    CPPUNIT_ASSERT(_submesh);
    Fields fields(*_submesh);

    const char* key = "field";
    const char* label = "displacement";
    fields.add(key, label);
    const size_t size = 1;
    CPPUNIT_ASSERT_EQUAL(size, fields._fields.size());

    PYLITH_METHOD_END;
} // testAdd


// ----------------------------------------------------------------------
// Test del().
void
pylith::topology::TestFieldsSubmesh::testDelete(void) {
    PYLITH_METHOD_BEGIN;

    CPPUNIT_ASSERT(_submesh);
    Fields fields(*_submesh);

    const char* keyA = "field A";
    const char* labelA = "displacement";
    fields.add(keyA, labelA);

    const char* keyB = "field B";
    const char* labelB = "velocity";
    fields.add(keyB, labelB);

    size_t size = 2;
    CPPUNIT_ASSERT_EQUAL(size, fields._fields.size());
    fields.del(keyA);
    size = 1;
    CPPUNIT_ASSERT_EQUAL(size, fields._fields.size());
    const Field& field = fields.get(keyB);
    CPPUNIT_ASSERT_EQUAL(std::string(labelB), std::string(field.getLabel()));

    PYLITH_METHOD_END;
} // testDelete


// ----------------------------------------------------------------------
// Test get().
void
pylith::topology::TestFieldsSubmesh::testGet(void) {
    PYLITH_METHOD_BEGIN;

    CPPUNIT_ASSERT(_submesh);
    Fields fields(*_submesh);

    const char* key = "field";
    const char* label = "displacement";
    fields.add(key, label);
    const Field& field = fields.get(key);
    CPPUNIT_ASSERT_EQUAL(std::string(label), std::string(field.getLabel()));

    PYLITH_METHOD_END;
} // testGet


// ----------------------------------------------------------------------
// Test get() const.
void
pylith::topology::TestFieldsSubmesh::testGetConst(void) {
    PYLITH_METHOD_BEGIN;

    CPPUNIT_ASSERT(_submesh);
    Fields fields(*_submesh);

    const char* key = "field";
    const char* label = "displacement";
    fields.add(key, label);

    const Fields* fieldsPtr = &fields;
    CPPUNIT_ASSERT(fieldsPtr);
    const Field& field = fieldsPtr->get(key);
    CPPUNIT_ASSERT_EQUAL(std::string(label), std::string(field.getLabel()));

    PYLITH_METHOD_END;
} // testGetConst


// ----------------------------------------------------------------------
// Test hasField().
void
pylith::topology::TestFieldsSubmesh::testHasField(void) {
    PYLITH_METHOD_BEGIN;

    CPPUNIT_ASSERT(_submesh);
    Fields fields(*_submesh);

    fields.add("field A", "velocity");

    CPPUNIT_ASSERT_EQUAL(true, fields.hasField("field A"));
    CPPUNIT_ASSERT_EQUAL(false, fields.hasField("field B"));
    CPPUNIT_ASSERT_EQUAL(false, fields.hasField("field C"));

    fields.add("field B", "displacement");

    CPPUNIT_ASSERT_EQUAL(true, fields.hasField("field A"));
    CPPUNIT_ASSERT_EQUAL(true, fields.hasField("field B"));
    CPPUNIT_ASSERT_EQUAL(false, fields.hasField("field C"));

    PYLITH_METHOD_END;
} // testHasField


// End of file
