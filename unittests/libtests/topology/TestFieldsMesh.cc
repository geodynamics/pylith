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
// Copyright (c) 2010-2012 University of California, Davis
//
// See COPYING for license information.
//
// ----------------------------------------------------------------------
//

#include <portinfo>

#include "TestFieldsMesh.hh" // Implementation of class methods

#include "pylith/topology/Mesh.hh" // USES Mesh
#include "pylith/topology/Field.hh" // USES Field
#include "pylith/topology/Fields.hh" // USES Fields

#include "pylith/meshio/MeshIOAscii.hh" // USES MeshIOAscii

// ----------------------------------------------------------------------
CPPUNIT_TEST_SUITE_REGISTRATION( pylith::topology::TestFieldsMesh );

// ----------------------------------------------------------------------
typedef pylith::topology::Fields<pylith::topology::Field<pylith::topology::Mesh> > FieldsMesh;
typedef pylith::topology::Mesh::RealSection RealSection;
typedef pylith::topology::Mesh::SieveMesh SieveMesh;

// ----------------------------------------------------------------------
void
pylith::topology::TestFieldsMesh::setUp(void)
{ // setUp
  _mesh = new Mesh;
  meshio::MeshIOAscii importer;
  importer.filename("data/tri3.mesh");
  importer.read(_mesh);
} // setUp

// ----------------------------------------------------------------------
void
pylith::topology::TestFieldsMesh::tearDown(void)
{ // tearDown
  delete _mesh; _mesh = 0;
} // tearDown

// ----------------------------------------------------------------------
// Test constructor.
void
pylith::topology::TestFieldsMesh::testConstructor(void)
{ // testConstructor
  CPPUNIT_ASSERT(0 != _mesh);
  FieldsMesh fields(*_mesh);
} // testConstructor
 
// ----------------------------------------------------------------------
// Test add().
void
pylith::topology::TestFieldsMesh::testAdd(void)
{ // testAdd
  CPPUNIT_ASSERT(0 != _mesh);
  FieldsMesh fields(*_mesh);
  
  const char* label = "field";
  fields.add(label, "displacement");
  const size_t size = 1;
  CPPUNIT_ASSERT_EQUAL(size, fields._fields.size());
} // testAdd

// ----------------------------------------------------------------------
// Test add(domain).
void
pylith::topology::TestFieldsMesh::testAddDomain(void)
{ // testAddDomain
  const int fiberDim = 3;

  CPPUNIT_ASSERT(0 != _mesh);
  FieldsMesh fields(*_mesh);
  
  const char* label = "field";
  fields.add(label, "velocity", Field<Mesh>::VERTICES_FIELD, fiberDim);
  const size_t size = 1;
  CPPUNIT_ASSERT_EQUAL(size, fields._fields.size());

  Field<Mesh>& field = fields.get(label);
  const ALE::Obj<RealSection>& section = field.section();
  CPPUNIT_ASSERT(!section.isNull());
  const ALE::Obj<SieveMesh>& sieveMesh = _mesh->sieveMesh();
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
pylith::topology::TestFieldsMesh::testDelete(void)
{ // testDelete
  CPPUNIT_ASSERT(0 != _mesh);
  FieldsMesh fields(*_mesh);

  const char* labelA = "field A";
  fields.add(labelA, "displacement");

  const char* labelB = "field B";
  fields.add(labelB, "velocity");

  size_t size = 2;
  CPPUNIT_ASSERT_EQUAL(size, fields._fields.size());
  fields.del(labelA);
  size = 1;
  CPPUNIT_ASSERT_EQUAL(size, fields._fields.size());
  const Field<Mesh>& field = fields.get(labelB);
  CPPUNIT_ASSERT_EQUAL(std::string("velocity"), std::string(field.label()));
} // testDelete

// ----------------------------------------------------------------------
// Test get().
void
pylith::topology::TestFieldsMesh::testGet(void)
{ // testGet
  CPPUNIT_ASSERT(0 != _mesh);
  FieldsMesh fields(*_mesh);

  const char* label = "field";
  fields.add(label, "velocity");
  const Field<Mesh>& field = fields.get(label);
} // testGet

// ----------------------------------------------------------------------
// Test get() const.
void
pylith::topology::TestFieldsMesh::testGetConst(void)
{ // testGetConst
  CPPUNIT_ASSERT(0 != _mesh);
  FieldsMesh fields(*_mesh);

  const char* label = "field";
  fields.add(label, "velocity");

  const FieldsMesh* fieldsPtr = &fields;
  CPPUNIT_ASSERT(0 != fieldsPtr);
  const Field<Mesh>& field = fieldsPtr->get(label);
} // testGetConst

// ----------------------------------------------------------------------
// Test hasField().
void
pylith::topology::TestFieldsMesh::testHasField(void)
{ // testHasField
  CPPUNIT_ASSERT(0 != _mesh);
  FieldsMesh fields(*_mesh);

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
pylith::topology::TestFieldsMesh::testCopyLayout(void)
{ // testCopyLayout
  const int fiberDim = 3;

  CPPUNIT_ASSERT(0 != _mesh);
  FieldsMesh fields(*_mesh);
  
  const char* labelA = "field A";
  fields.add(labelA, "displacement", Field<Mesh>::VERTICES_FIELD, fiberDim);

  const char* labelB = "field B";
  fields.add(labelB, "velocity");
  Field<Mesh>& fieldA = fields.get(labelA);
  fieldA.allocate();

  fields.copyLayout(labelA);

  const size_t size = 2;
  CPPUNIT_ASSERT_EQUAL(size, fields._fields.size());
  const Field<Mesh>& field = fields.get(labelB);
  const ALE::Obj<RealSection>& section = field.section();
  CPPUNIT_ASSERT(!section.isNull());
  const ALE::Obj<SieveMesh>& sieveMesh = _mesh->sieveMesh();
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
