// -*- C++ -*-
//
// ----------------------------------------------------------------------
//
//                           Brad T. Aagaard
//                        U.S. Geological Survey
//
// {LicenseText}
//
// ----------------------------------------------------------------------
//

#include <portinfo>

#include "TestFieldsMesh.hh" // Implementation of class methods

#include "pylith/topology/Fields.hh" // USES Fields

#include "pylith/topology/Field.hh" // USES Field
#include "pylith/topology/Mesh.hh" // USES Mesh
#include "pylith/meshio/MeshIOAscii.hh" // USES MeshIOAscii

typedef pylith::topology::Fields<pylith::topology::Field,
				 pylith::topology::Mesh> FieldsMesh;

// ----------------------------------------------------------------------
CPPUNIT_TEST_SUITE_REGISTRATION( pylith::topology::TestFieldsMesh );

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
  fields.add(label);
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
  fields.add(label, Field::VERTICES_FIELD, fiberDim);
  const size_t size = 1;
  CPPUNIT_ASSERT_EQUAL(size, fields._fields.size());

  Field& field = fields.get(label);
  const ALE::Obj<MeshRealSection>& section = field.section();
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
  fields.add(labelA);

  const char* labelB = "field B";
  fields.add(labelB);

  size_t size = 2;
  CPPUNIT_ASSERT_EQUAL(size, fields._fields.size());
  fields.del(labelA);
  size = 1;
  CPPUNIT_ASSERT_EQUAL(size, fields._fields.size());
  const Field& field = fields.get(labelB);
} // testDelete

// ----------------------------------------------------------------------
// Test get().
void
pylith::topology::TestFieldsMesh::testGet(void)
{ // testGet
  CPPUNIT_ASSERT(0 != _mesh);
  FieldsMesh fields(*_mesh);

  const char* label = "field";
  fields.add(label);
  const Field& field = fields.get(label);
} // testGet

// ----------------------------------------------------------------------
// Test get() const.
void
pylith::topology::TestFieldsMesh::testGetConst(void)
{ // testGetConst
  CPPUNIT_ASSERT(0 != _mesh);
  FieldsMesh fields(*_mesh);

  const char* label = "field";
  fields.add(label);

  const FieldsMesh* fieldsPtr = &fields;
  CPPUNIT_ASSERT(0 != fieldsPtr);
  const Field& field = fieldsPtr->get(label);
} // testGetConst

// ----------------------------------------------------------------------
// Test copyLayout().
void
pylith::topology::TestFieldsMesh::testCopyLayout(void)
{ // testCopyLayout
  const int fiberDim = 3;

  CPPUNIT_ASSERT(0 != _mesh);
  FieldsMesh fields(*_mesh);
  
  const char* labelA = "field A";
  fields.add(labelA, Field::VERTICES_FIELD, fiberDim);

  const char* labelB = "field B";
  fields.add(labelB);
  Field& fieldA = fields.get(labelA);
  fieldA.allocate();

  fields.copyLayout(labelA);

  const size_t size = 2;
  CPPUNIT_ASSERT_EQUAL(size, fields._fields.size());
  const Field& field = fields.get(labelB);
  const ALE::Obj<MeshRealSection>& section = field.section();
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
