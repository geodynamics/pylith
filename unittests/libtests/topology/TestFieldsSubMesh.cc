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
// Copyright (c) 2010-2013 University of California, Davis
//
// See COPYING for license information.
//
// ----------------------------------------------------------------------
//

#include <portinfo>

#include "TestFieldsSubMesh.hh" // Implementation of class methods

#include "pylith/topology/SubMesh.hh" // USES SubMesh
#include "pylith/topology/Field.hh" // USES Field
#include "pylith/topology/Fields.hh" // USES Fields
#include "pylith/topology/Stratum.hh" // USES Stratum
#include "pylith/topology/VisitorMesh.hh" // USES VecVisitorMesh

#include "pylith/meshio/MeshIOAscii.hh" // USES MeshIOAscii

// ----------------------------------------------------------------------
CPPUNIT_TEST_SUITE_REGISTRATION( pylith::topology::TestFieldsSubMesh );

// ----------------------------------------------------------------------
typedef pylith::topology::Fields<pylith::topology::Field<pylith::topology::SubMesh> > FieldsSubMesh;

// ----------------------------------------------------------------------
void
pylith::topology::TestFieldsSubMesh::setUp(void)
{ // setUp
  PYLITH_METHOD_BEGIN;

  _mesh = new Mesh;
  meshio::MeshIOAscii importer;
  importer.filename("data/tri3.mesh");
  importer.read(_mesh);

  _submesh = new SubMesh(*_mesh, "bc");

  PYLITH_METHOD_END;
} // setUp

// ----------------------------------------------------------------------
void
pylith::topology::TestFieldsSubMesh::tearDown(void)
{ // tearDown
  PYLITH_METHOD_BEGIN;

  delete _mesh; _mesh = 0;
  delete _submesh; _submesh = 0;

  PYLITH_METHOD_END;
} // tearDown

// ----------------------------------------------------------------------
// Test constructor.
void
pylith::topology::TestFieldsSubMesh::testConstructor(void)
{ // testConstructor
  PYLITH_METHOD_BEGIN;

  CPPUNIT_ASSERT(_submesh);
  FieldsSubMesh fields(*_submesh);

  PYLITH_METHOD_END;
} // testConstructor
 
// ----------------------------------------------------------------------
// Test add().
void
pylith::topology::TestFieldsSubMesh::testAdd(void)
{ // testAdd
  PYLITH_METHOD_BEGIN;

  CPPUNIT_ASSERT(_submesh);
  FieldsSubMesh fields(*_submesh);
  
  const char* label = "field";
  fields.add(label, "displacement");
  const size_t size = 1;
  CPPUNIT_ASSERT_EQUAL(size, fields._fields.size());

  PYLITH_METHOD_END;
} // testAdd

// ----------------------------------------------------------------------
// Test add(domain).
void
pylith::topology::TestFieldsSubMesh::testAddDomain(void)
{ // testAddDomain
  PYLITH_METHOD_BEGIN;

  const int fiberDim = 3;

  CPPUNIT_ASSERT(_submesh);
  FieldsSubMesh fields(*_submesh);

  const char* label = "field";
  fields.add(label, "velocity", Field<SubMesh>::VERTICES_FIELD, fiberDim);
  const size_t size = 1;
  CPPUNIT_ASSERT_EQUAL(size, fields._fields.size());

  Field<SubMesh>& field = fields.get(label);
  field.allocate();

  PetscDM dmMesh = _submesh->dmMesh();CPPUNIT_ASSERT(dmMesh);
  Stratum depthStratum(dmMesh, Stratum::DEPTH, 0);
  const PetscInt vStart = depthStratum.begin();
  const PetscInt vEnd = depthStratum.end();

  VecVisitorMesh fieldVisitor(field);
  for(PetscInt v = vStart; v < vEnd; ++v) {
    CPPUNIT_ASSERT_EQUAL(fiberDim, fieldVisitor.sectionDof(v));
  } // for

  PYLITH_METHOD_END;
} // testAddDomain

// ----------------------------------------------------------------------
// Test del().
void
pylith::topology::TestFieldsSubMesh::testDelete(void)
{ // testDelete
  PYLITH_METHOD_BEGIN;

  CPPUNIT_ASSERT(_submesh);
  FieldsSubMesh fields(*_submesh);

  const char* labelA = "field A";
  fields.add(labelA, "displacement");

  const char* labelB = "field B";
  fields.add(labelB, "velocity");

  size_t size = 2;
  CPPUNIT_ASSERT_EQUAL(size, fields._fields.size());
  fields.del(labelA);
  size = 1;
  CPPUNIT_ASSERT_EQUAL(size, fields._fields.size());
  const Field<SubMesh>& field = fields.get(labelB);

  PYLITH_METHOD_END;
} // testDelete

// ----------------------------------------------------------------------
// Test get().
void
pylith::topology::TestFieldsSubMesh::testGet(void)
{ // testGet
  PYLITH_METHOD_BEGIN;

  CPPUNIT_ASSERT(_submesh);
  FieldsSubMesh fields(*_submesh);

  const char* label = "field";
  fields.add(label, "displacement");
  const Field<SubMesh>& field = fields.get(label);

  PYLITH_METHOD_END;
} // testGet

// ----------------------------------------------------------------------
// Test get() const.
void
pylith::topology::TestFieldsSubMesh::testGetConst(void)
{ // testGetConst
  PYLITH_METHOD_BEGIN;

  CPPUNIT_ASSERT(_submesh);
  FieldsSubMesh fields(*_submesh);

  const char* label = "field";
  fields.add(label, "displacement");

  const FieldsSubMesh* fieldsPtr = &fields;
  CPPUNIT_ASSERT(fieldsPtr);
  const Field<SubMesh>& field = fieldsPtr->get(label);

  PYLITH_METHOD_END;
} // testGetConst

// ----------------------------------------------------------------------
// Test hasField().
void
pylith::topology::TestFieldsSubMesh::testHasField(void)
{ // testHasField
  PYLITH_METHOD_BEGIN;

  CPPUNIT_ASSERT(_submesh);
  FieldsSubMesh fields(*_submesh);

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

// ----------------------------------------------------------------------
// Test copyLayout().
void
pylith::topology::TestFieldsSubMesh::testCopyLayout(void)
{ // testCopyLayout
  PYLITH_METHOD_BEGIN;

  const int fiberDim = 3;

  CPPUNIT_ASSERT(_submesh);
  FieldsSubMesh fields(*_submesh);
  
  const char* labelA = "field A";
  fields.add(labelA, "velocity", Field<SubMesh>::VERTICES_FIELD, fiberDim);

  const char* labelB = "field B";
  fields.add(labelB, "displacement");
  Field<SubMesh>& fieldA = fields.get(labelA);
  fieldA.allocate();

  fields.copyLayout(labelA);

  const size_t size = 2;
  CPPUNIT_ASSERT_EQUAL(size, fields._fields.size());
  const Field<SubMesh>& field = fields.get(labelB);

  PetscDM dmMesh = _submesh->dmMesh();CPPUNIT_ASSERT(dmMesh);
  Stratum depthStratum(dmMesh, Stratum::DEPTH, 0);
  const PetscInt vStart = depthStratum.begin();
  const PetscInt vEnd = depthStratum.end();

  VecVisitorMesh fieldVisitor(field);
  for(PetscInt v = vStart; v < vEnd; ++v) {
    CPPUNIT_ASSERT_EQUAL(fiberDim, fieldVisitor.sectionDof(v));
  } // for

  PYLITH_METHOD_END;
} // testCopyLayout


// End of file 
