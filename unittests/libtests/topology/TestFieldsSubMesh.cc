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

#include "pylith/meshio/MeshIOAscii.hh" // USES MeshIOAscii

// ----------------------------------------------------------------------
CPPUNIT_TEST_SUITE_REGISTRATION( pylith::topology::TestFieldsSubMesh );

// ----------------------------------------------------------------------
typedef pylith::topology::Fields<pylith::topology::Field<pylith::topology::SubMesh> > FieldsSubMesh;
typedef pylith::topology::SubMesh::RealSection RealSection;
typedef pylith::topology::SubMesh::SieveMesh SieveMesh;

// ----------------------------------------------------------------------
void
pylith::topology::TestFieldsSubMesh::setUp(void)
{ // setUp
  _mesh = new Mesh;
  meshio::MeshIOAscii importer;
  importer.filename("data/tri3.mesh");
  importer.read(_mesh);

  _submesh = new SubMesh(*_mesh, "bc");
} // setUp

// ----------------------------------------------------------------------
void
pylith::topology::TestFieldsSubMesh::tearDown(void)
{ // tearDown
  delete _mesh; _mesh = 0;
  delete _submesh; _submesh = 0;
} // tearDown

// ----------------------------------------------------------------------
// Test constructor.
void
pylith::topology::TestFieldsSubMesh::testConstructor(void)
{ // testConstructor
  CPPUNIT_ASSERT(0 != _submesh);
  FieldsSubMesh fields(*_submesh);
} // testConstructor
 
// ----------------------------------------------------------------------
// Test add().
void
pylith::topology::TestFieldsSubMesh::testAdd(void)
{ // testAdd
  CPPUNIT_ASSERT(0 != _submesh);
  FieldsSubMesh fields(*_submesh);
  
  const char* label = "field";
  fields.add(label, "displacement");
  const size_t size = 1;
  CPPUNIT_ASSERT_EQUAL(size, fields._fields.size());
} // testAdd

// ----------------------------------------------------------------------
// Test add(domain).
void
pylith::topology::TestFieldsSubMesh::testAddDomain(void)
{ // testAddDomain
  const int fiberDim = 3;

  CPPUNIT_ASSERT(0 != _submesh);
  FieldsSubMesh fields(*_submesh);
  
  const char* label = "field";
  fields.add(label, "velocity", Field<SubMesh>::VERTICES_FIELD, fiberDim);
  const size_t size = 1;
  CPPUNIT_ASSERT_EQUAL(size, fields._fields.size());

  Field<SubMesh>& field = fields.get(label);
  const ALE::Obj<RealSection>& section = field.section();
  CPPUNIT_ASSERT(!section.isNull());
  const ALE::Obj<SieveMesh>& sieveMesh = _submesh->sieveMesh();
  CPPUNIT_ASSERT(!sieveMesh.isNull());
  const ALE::Obj<SieveMesh::label_sequence>& vertices = 
    sieveMesh->depthStratum(0);
  CPPUNIT_ASSERT(!vertices.isNull());
  field.allocate();
  for (SieveMesh::label_sequence::iterator v_iter=vertices->begin();
       v_iter != vertices->end();
       ++v_iter)
    CPPUNIT_ASSERT_EQUAL(fiberDim, section->getFiberDimension(*v_iter));
} // testAddDomain

// ----------------------------------------------------------------------
// Test del().
void
pylith::topology::TestFieldsSubMesh::testDelete(void)
{ // testDelete
  CPPUNIT_ASSERT(0 != _submesh);
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
} // testDelete

// ----------------------------------------------------------------------
// Test get().
void
pylith::topology::TestFieldsSubMesh::testGet(void)
{ // testGet
  CPPUNIT_ASSERT(0 != _submesh);
  FieldsSubMesh fields(*_submesh);

  const char* label = "field";
  fields.add(label, "displacement");
  const Field<SubMesh>& field = fields.get(label);
} // testGet

// ----------------------------------------------------------------------
// Test get() const.
void
pylith::topology::TestFieldsSubMesh::testGetConst(void)
{ // testGetConst
  CPPUNIT_ASSERT(0 != _submesh);
  FieldsSubMesh fields(*_submesh);

  const char* label = "field";
  fields.add(label, "displacement");

  const FieldsSubMesh* fieldsPtr = &fields;
  CPPUNIT_ASSERT(0 != fieldsPtr);
  const Field<SubMesh>& field = fieldsPtr->get(label);
} // testGetConst

// ----------------------------------------------------------------------
// Test hasField().
void
pylith::topology::TestFieldsSubMesh::testHasField(void)
{ // testHasField
  CPPUNIT_ASSERT(0 != _submesh);
  FieldsSubMesh fields(*_submesh);

  fields.add("field A", "velocity");
  
  CPPUNIT_ASSERT_EQUAL(true, fields.hasField("field A"));
  CPPUNIT_ASSERT_EQUAL(false, fields.hasField("field B"));
  CPPUNIT_ASSERT_EQUAL(false, fields.hasField("field C"));

  fields.add("field B", "displacement");

  CPPUNIT_ASSERT_EQUAL(true, fields.hasField("field A"));
  CPPUNIT_ASSERT_EQUAL(true, fields.hasField("field B"));
  CPPUNIT_ASSERT_EQUAL(false, fields.hasField("field C"));

} // testHasField

// ----------------------------------------------------------------------
// Test copyLayout().
void
pylith::topology::TestFieldsSubMesh::testCopyLayout(void)
{ // testCopyLayout
  const int fiberDim = 3;

  CPPUNIT_ASSERT(0 != _submesh);
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
  const ALE::Obj<RealSection>& section = field.section();
  CPPUNIT_ASSERT(!section.isNull());
  const ALE::Obj<SieveMesh>& sieveMesh = _submesh->sieveMesh();
  CPPUNIT_ASSERT(!sieveMesh.isNull());
  const ALE::Obj<SieveMesh::label_sequence>& vertices
    = sieveMesh->depthStratum(0);
  CPPUNIT_ASSERT(!vertices.isNull());
  for (SieveMesh::label_sequence::iterator v_iter=vertices->begin();
       v_iter != vertices->end();
       ++v_iter)
    CPPUNIT_ASSERT_EQUAL(fiberDim, section->getFiberDimension(*v_iter));
} // testCopyLayout


// End of file 
