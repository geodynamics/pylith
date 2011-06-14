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
// Copyright (c) 2010-2011 University of California, Davis
//
// See COPYING for license information.
//
// ----------------------------------------------------------------------
//

#include <portinfo>

#include "TestFieldsNewMesh.hh" // Implementation of class methods

#include "pylith/topology/Mesh.hh" // USES Mesh
#include "pylith/topology/Field.hh" // USES Field
#include "pylith/topology/FieldsNew.hh" // USES FieldsNew

#include "pylith/meshio/MeshIOAscii.hh" // USES MeshIOAscii

// ----------------------------------------------------------------------
CPPUNIT_TEST_SUITE_REGISTRATION( pylith::topology::TestFieldsNewMesh );

// ----------------------------------------------------------------------
typedef pylith::topology::FieldsNew<pylith::topology::Mesh> FieldsNewMesh;
typedef pylith::topology::Mesh::SieveMesh SieveMesh;
typedef pylith::topology::Mesh::RealUniformSection RealUniformSection;

// ----------------------------------------------------------------------
void
pylith::topology::TestFieldsNewMesh::setUp(void)
{ // setUp
  _mesh = new Mesh;
  meshio::MeshIOAscii importer;
  importer.filename("data/tri3.mesh");
  importer.read(_mesh);
} // setUp

// ----------------------------------------------------------------------
void
pylith::topology::TestFieldsNewMesh::tearDown(void)
{ // tearDown
  delete _mesh; _mesh = 0;
} // tearDown

// ----------------------------------------------------------------------
// Test constructor.
void
pylith::topology::TestFieldsNewMesh::testConstructor(void)
{ // testConstructor
  CPPUNIT_ASSERT(0 != _mesh);
  FieldsNewMesh fields(*_mesh);
} // testConstructor
 
// ----------------------------------------------------------------------
// Test hasField().
void
pylith::topology::TestFieldsNewMesh::testHasField(void)
{ // testHasField
  CPPUNIT_ASSERT(0 != _mesh);
  FieldsNewMesh fields(*_mesh);

  fields.add("field A", "velocity", 3, FieldBase::VECTOR);
  
  CPPUNIT_ASSERT_EQUAL(true, fields.hasField("field A"));
  CPPUNIT_ASSERT_EQUAL(false, fields.hasField("field B"));
  CPPUNIT_ASSERT_EQUAL(false, fields.hasField("field C"));

  fields.add("field B", "displacement", 3, FieldBase::VECTOR);

  CPPUNIT_ASSERT_EQUAL(true, fields.hasField("field A"));
  CPPUNIT_ASSERT_EQUAL(true, fields.hasField("field B"));
  CPPUNIT_ASSERT_EQUAL(false, fields.hasField("field C"));

} // testHasField

// ----------------------------------------------------------------------
// Test add().
void
pylith::topology::TestFieldsNewMesh::testAdd(void)
{ // testAdd
  CPPUNIT_ASSERT(0 != _mesh);
  FieldsNewMesh fields(*_mesh);
  
  fields.add("field A", "velocity", 3, FieldBase::VECTOR);
  fields.add("field B", "displacement", 4, FieldBase::OTHER, 2.0, true);

  const size_t size = 2;
  CPPUNIT_ASSERT_EQUAL(size, fields._fields.size());
} // testAdd

// ----------------------------------------------------------------------
// Test allocate(sequence).
void
pylith::topology::TestFieldsNewMesh::testAllocateSequence(void)
{ // testAllocateSequence
  const int fiberDim = 7;

  CPPUNIT_ASSERT(0 != _mesh);
  FieldsNewMesh fields(*_mesh);
  
  fields.add("field A", "velocity", 3, FieldBase::VECTOR);
  fields.add("field B", "other", 4, FieldBase::OTHER);

  const ALE::Obj<SieveMesh>& sieveMesh = _mesh->sieveMesh();
  CPPUNIT_ASSERT(!sieveMesh.isNull());
  const ALE::Obj<SieveMesh::label_sequence>& vertices = 
    sieveMesh->depthStratum(0);
  CPPUNIT_ASSERT(!vertices.isNull());

  fields.allocate(vertices);

  const size_t size = 2;
  CPPUNIT_ASSERT_EQUAL(size, fields._fields.size());

  const ALE::Obj<RealUniformSection>& section = fields.section();
  CPPUNIT_ASSERT(!section.isNull());
  for (SieveMesh::label_sequence::iterator v_iter=vertices->begin();
       v_iter != vertices->end();
       ++v_iter)
    CPPUNIT_ASSERT_EQUAL(fiberDim, section->getFiberDimension(*v_iter));
} // testAllocateSequence

// ----------------------------------------------------------------------
// Test allocate(array).
void
pylith::topology::TestFieldsNewMesh::testAllocateArray(void)
{ // testAllocateSequence
  const int fiberDim = 7;
  const int nptsIn = 3;
  const int ptsIn[nptsIn] = {
    1, 3, 4,
  };
  const int nptsOut = 1;
  const int ptsOut[nptsOut] = {
    2,
  };

  CPPUNIT_ASSERT(0 != _mesh);
  FieldsNewMesh fields(*_mesh);
  
  fields.add("field A", "velocity", 3, FieldBase::VECTOR);
  fields.add("field B", "other", 4, FieldBase::OTHER);

  int_array verticesIn(nptsIn);
  for (int i=0; i < nptsIn; ++i)
    verticesIn = ptsIn[i];

  int_array verticesOut(nptsOut);
  for (int i=0; i < nptsOut; ++i)
    verticesOut = ptsOut[i];

  fields.allocate(verticesIn);

  const size_t size = 2;
  CPPUNIT_ASSERT_EQUAL(size, fields._fields.size());

  const ALE::Obj<RealUniformSection>& section = fields.section();
  CPPUNIT_ASSERT(!section.isNull());
  for (int i=0; i < nptsIn; ++i)
    CPPUNIT_ASSERT_EQUAL(fiberDim, section->getFiberDimension(verticesIn[i]));
  for (int i=0; i < nptsOut; ++i)
    CPPUNIT_ASSERT_EQUAL(0, section->getFiberDimension(verticesOut[i]));
} // testAllocateArray

// ----------------------------------------------------------------------
// Test allocate(domain).
void
pylith::topology::TestFieldsNewMesh::testAllocateDomain(void)
{ // testAllocateDomain
  const int fiberDim = 7;

  CPPUNIT_ASSERT(0 != _mesh);
  FieldsNewMesh fields(*_mesh);
  
  fields.add("field A", "velocity", 3, FieldBase::VECTOR);
  fields.add("field B", "other", 4, FieldBase::OTHER);
  fields.allocate(Field<Mesh>::VERTICES_FIELD);

  const size_t size = 2;
  CPPUNIT_ASSERT_EQUAL(size, fields._fields.size());

  const ALE::Obj<RealUniformSection>& section = fields.section();
  CPPUNIT_ASSERT(!section.isNull());
  const ALE::Obj<SieveMesh>& sieveMesh = _mesh->sieveMesh();
  CPPUNIT_ASSERT(!sieveMesh.isNull());
  const ALE::Obj<SieveMesh::label_sequence>& vertices = 
    sieveMesh->depthStratum(0);
  CPPUNIT_ASSERT(!vertices.isNull());
  for (SieveMesh::label_sequence::iterator v_iter=vertices->begin();
       v_iter != vertices->end();
       ++v_iter)
    CPPUNIT_ASSERT_EQUAL(fiberDim, section->getFiberDimension(*v_iter));
} // testAllocateDomain

// ----------------------------------------------------------------------
// Test get().
void
pylith::topology::TestFieldsNewMesh::testGet(void)
{ // testGet
  CPPUNIT_ASSERT(0 != _mesh);
  FieldsNewMesh fields(*_mesh);

  fields.add("field A", "velocity", 3, FieldBase::VECTOR);
  fields.add("field B", "displacement", 4, FieldBase::OTHER, 2.0, true);
  fields.allocate(FieldBase::VERTICES_FIELD);

  const ALE::Obj<Mesh::SieveMesh>& sieveMesh = _mesh->sieveMesh();
  CPPUNIT_ASSERT(!sieveMesh.isNull());
  const ALE::Obj<Mesh::SieveMesh::label_sequence>& vertices = 
    sieveMesh->depthStratum(0);
  CPPUNIT_ASSERT(!vertices.isNull());
  const Mesh::SieveMesh::label_sequence::iterator verticesBegin =
    vertices->begin();
  const Mesh::SieveMesh::label_sequence::iterator verticesEnd =
    vertices->end();

  // Check field A
  Field<Mesh>& fieldA = fields.get("field A");
  CPPUNIT_ASSERT_EQUAL(std::string("velocity"), std::string(fieldA.label()));
  CPPUNIT_ASSERT_EQUAL(FieldBase::VECTOR,
		       fieldA.vectorFieldType());
  CPPUNIT_ASSERT_EQUAL(1.0, fieldA.scale());
  CPPUNIT_ASSERT_EQUAL(false, fieldA.addDimensionOkay());

  const ALE::Obj<Mesh::RealSection>& sectionA = fieldA.section();
  CPPUNIT_ASSERT(!sectionA.isNull());
  for(Mesh::SieveMesh::label_sequence::iterator v_iter = verticesBegin;
      v_iter != verticesEnd;
      ++v_iter) {
    const int fiberDim = sectionA->getFiberDimension(*v_iter);
    CPPUNIT_ASSERT_EQUAL(3, fiberDim);
    const double* values = sectionA->restrictPoint(*v_iter);
    CPPUNIT_ASSERT(values);
  } // for


  // Check field B
  Field<Mesh>& fieldB = fields.get("field B");
  CPPUNIT_ASSERT_EQUAL(std::string("displacement"), 
		       std::string(fieldB.label()));
  CPPUNIT_ASSERT_EQUAL(FieldBase::OTHER,
		       fieldB.vectorFieldType());
  CPPUNIT_ASSERT_EQUAL(2.0, fieldB.scale());
  CPPUNIT_ASSERT_EQUAL(true, fieldB.addDimensionOkay());

  const ALE::Obj<Mesh::RealSection>& sectionB = fieldB.section();
  CPPUNIT_ASSERT(!sectionB.isNull());
  for(Mesh::SieveMesh::label_sequence::iterator v_iter = verticesBegin;
      v_iter != verticesEnd;
      ++v_iter) {
    const int fiberDim = sectionB->getFiberDimension(*v_iter);
    CPPUNIT_ASSERT_EQUAL(4, fiberDim);
    const double* values = sectionB->restrictPoint(*v_iter);
    CPPUNIT_ASSERT(values);
  } // for

  // Make sure we can clone field B
  Field<Mesh> fieldC(*_mesh);
  fieldC.cloneSection(fieldB);
  fieldC.copy(fieldB);
} // testGet

// ----------------------------------------------------------------------
// Test get() const.
void
pylith::topology::TestFieldsNewMesh::testGetConst(void)
{ // testGetConst
  CPPUNIT_ASSERT(0 != _mesh);
  FieldsNewMesh fields(*_mesh);

  fields.add("field A", "velocity", 3, FieldBase::VECTOR);
  fields.add("field B", "displacement", 4, FieldBase::OTHER, 2.0, true);
  fields.allocate(FieldBase::VERTICES_FIELD);

  const ALE::Obj<Mesh::SieveMesh>& sieveMesh = _mesh->sieveMesh();
  CPPUNIT_ASSERT(!sieveMesh.isNull());
  const ALE::Obj<Mesh::SieveMesh::label_sequence>& vertices = 
    sieveMesh->depthStratum(0);
  CPPUNIT_ASSERT(!vertices.isNull());
  const Mesh::SieveMesh::label_sequence::iterator verticesBegin =
    vertices->begin();
  const Mesh::SieveMesh::label_sequence::iterator verticesEnd =
    vertices->end();

  // Check field A
  const Field<Mesh>& fieldA = fields.get("field A");
  CPPUNIT_ASSERT_EQUAL(std::string("velocity"), std::string(fieldA.label()));
  CPPUNIT_ASSERT_EQUAL(FieldBase::VECTOR,
		       fieldA.vectorFieldType());
  CPPUNIT_ASSERT_EQUAL(1.0, fieldA.scale());
  CPPUNIT_ASSERT_EQUAL(false, fieldA.addDimensionOkay());

  const ALE::Obj<Mesh::RealSection>& sectionA = fieldA.section();
  CPPUNIT_ASSERT(!sectionA.isNull());
  for(Mesh::SieveMesh::label_sequence::iterator v_iter = verticesBegin;
      v_iter != verticesEnd;
      ++v_iter) {
    const int fiberDim = sectionA->getFiberDimension(*v_iter);
    CPPUNIT_ASSERT_EQUAL(3, fiberDim);
  } // for


  // Check field B
  const Field<Mesh>& fieldB = fields.get("field B");
  CPPUNIT_ASSERT_EQUAL(std::string("displacement"), 
		       std::string(fieldB.label()));
  CPPUNIT_ASSERT_EQUAL(FieldBase::OTHER,
		       fieldB.vectorFieldType());
  CPPUNIT_ASSERT_EQUAL(2.0, fieldB.scale());
  CPPUNIT_ASSERT_EQUAL(true, fieldB.addDimensionOkay());

  const ALE::Obj<Mesh::RealSection>& sectionB = fieldB.section();
  CPPUNIT_ASSERT(!sectionB.isNull());
  for(Mesh::SieveMesh::label_sequence::iterator v_iter = verticesBegin;
      v_iter != verticesEnd;
      ++v_iter) {
    const int fiberDim = sectionB->getFiberDimension(*v_iter);
    CPPUNIT_ASSERT_EQUAL(4, fiberDim);
  } // for
} // testGetConst

// ----------------------------------------------------------------------
// Test mesh().
void
pylith::topology::TestFieldsNewMesh::testMesh(void)
{ // testMesh
  CPPUNIT_ASSERT(0 != _mesh);
  FieldsNewMesh fields(*_mesh);

  const Mesh& mesh = fields.mesh();
} // testMesh

// ----------------------------------------------------------------------
// Test fiberDim().
void
pylith::topology::TestFieldsNewMesh::testFiberDim(void)
{ // testFiberDim
  CPPUNIT_ASSERT(0 != _mesh);
  FieldsNewMesh fields(*_mesh);

  fields.add("field A", "velocity", 3, FieldBase::VECTOR);
  fields.add("field B", "displacement", 4, FieldBase::OTHER, 2.0, true);

  CPPUNIT_ASSERT_EQUAL(7, fields.fiberDim());
} // testFiberDim

// ----------------------------------------------------------------------
// Test sectionIndex().
void
pylith::topology::TestFieldsNewMesh::testSectionIndex(void)
{ // testSectionIndex
  CPPUNIT_ASSERT(0 != _mesh);
  FieldsNewMesh fields(*_mesh);

  fields.add("field A", "velocity", 3, FieldBase::VECTOR);
  fields.add("field B", "displacement", 4, FieldBase::OTHER, 2.0, true);

  CPPUNIT_ASSERT_EQUAL(0, fields.sectionIndex("field A"));
  CPPUNIT_ASSERT_EQUAL(3, fields.sectionIndex("field B"));
} // testSectionIndex

// ----------------------------------------------------------------------
// Test sectionFiberDim().
void
pylith::topology::TestFieldsNewMesh::testSectionFiberDim(void)
{ // testSectionFiberDim
  CPPUNIT_ASSERT(0 != _mesh);
  FieldsNewMesh fields(*_mesh);

  fields.add("field A", "velocity", 3, FieldBase::VECTOR);
  fields.add("field B", "displacement", 4, FieldBase::OTHER, 2.0, true);

  CPPUNIT_ASSERT_EQUAL(3, fields.sectionFiberDim("field A"));
  CPPUNIT_ASSERT_EQUAL(4, fields.sectionFiberDim("field B"));
} // testSectionFiberDim

// ----------------------------------------------------------------------
// Test complete().
void
pylith::topology::TestFieldsNewMesh::testComplete(void)
{ // testComplete
  CPPUNIT_ASSERT(0 != _mesh);
  FieldsNewMesh fields(*_mesh);

  fields.add("field A", "velocity", 3, FieldBase::VECTOR);
  fields.add("field B", "displacement", 4, FieldBase::OTHER, 2.0, true);
  fields.allocate(FieldBase::VERTICES_FIELD);

  fields.complete();
} // testComplete

// ----------------------------------------------------------------------
// Test fieldNames() const.
void
pylith::topology::TestFieldsNewMesh::testFieldNames(void)
{ // testFieldNames
  const int numFieldsE = 2;
  const char* namesE[2] = {
    "field A",
    "field B"
  };

  CPPUNIT_ASSERT(0 != _mesh);
  FieldsNewMesh fields(*_mesh);

  fields.add("field A", "velocity", 3, FieldBase::VECTOR);
  fields.add("field B", "displacement", 4, FieldBase::OTHER, 2.0, true);

  int numFields = 0;
  char** names = 0;
  fields.fieldNames(&numFields, &names);
  
  CPPUNIT_ASSERT_EQUAL(numFieldsE, numFields);
  
  for (int i=0; i < numFields; ++i)
    CPPUNIT_ASSERT_EQUAL(std::string(namesE[i]), std::string(names[i]));

  for (int i=0; i < numFields; ++i) {
    delete[] names[i]; names[i] = 0;
  } // for
  delete[] names; names = 0;
} // testFieldNames

// ----------------------------------------------------------------------
// Test view().
void
pylith::topology::TestFieldsNewMesh::testView(void)
{ // testView
  CPPUNIT_ASSERT(0 != _mesh);
  FieldsNewMesh fields(*_mesh);

  fields.add("field A", "velocity", 3, FieldBase::VECTOR);
  fields.add("field B", "displacement", 4, FieldBase::OTHER, 2.0, true);
  fields.allocate(FieldBase::VERTICES_FIELD);

  fields.view("TEST VIEW");
} // testView


// End of file 
