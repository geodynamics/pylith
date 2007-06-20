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

#include "TestFieldsManager.hh" // Implementation of class methods

#include "pylith/topology/FieldsManager.hh" // USES FieldsManager

#include "pylith/meshio/MeshIOAscii.hh" // USES MeshIOAscii
#include "pylith/utils/sievetypes.hh" // USES PETSc Mesh

// ----------------------------------------------------------------------
CPPUNIT_TEST_SUITE_REGISTRATION( pylith::topology::TestFieldsManager );

// ----------------------------------------------------------------------
// Test constructor.
void
pylith::topology::TestFieldsManager::testConstructor(void)
{ // testConstructor
  ALE::Obj<Mesh> mesh;
  _initialize(&mesh);
  FieldsManager manager(mesh);
} // testConstructor

// ----------------------------------------------------------------------
// Test addReal().
void
pylith::topology::TestFieldsManager::testAddReal(void)
{ // testAddReal
  ALE::Obj<Mesh> mesh;
  _initialize(&mesh);
  FieldsManager manager(mesh);

  const char* label = "field";
  manager.addReal(label);
  const size_t size = 1;
  CPPUNIT_ASSERT_EQUAL(size, manager._real.size());
} // testAddReal

// ----------------------------------------------------------------------
// Test getReal().
void
pylith::topology::TestFieldsManager::testGetReal(void)
{ // testGetReal
  ALE::Obj<Mesh> mesh;
  _initialize(&mesh);
  FieldsManager manager(mesh);

  const char* label = "field";
  manager.addReal(label);
  const ALE::Obj<real_section_type>& field = manager.getReal(label);
  CPPUNIT_ASSERT(!field.isNull());
} // testGetReal

// ----------------------------------------------------------------------
// Test delReal().
void
pylith::topology::TestFieldsManager::testDelReal(void)
{ // testDelReal
  ALE::Obj<Mesh> mesh;
  _initialize(&mesh);
  FieldsManager manager(mesh);

  const char* labelA = "field A";
  manager.addReal(labelA);

  const char* labelB = "field B";
  manager.addReal(labelB);

  size_t size = 2;
  CPPUNIT_ASSERT_EQUAL(size, manager._real.size());
  manager.delReal(labelA);
  size = 1;
  CPPUNIT_ASSERT_EQUAL(size, manager._real.size());
  const ALE::Obj<real_section_type>& field = manager.getReal(labelB);
  CPPUNIT_ASSERT(!field.isNull());
} // testDelReal

// ----------------------------------------------------------------------
// Test setFiberDimension().
void
pylith::topology::TestFieldsManager::testSetFiberDimension(void)
{ // testSetFiberDimension
  ALE::Obj<Mesh> mesh;
  _initialize(&mesh);
  FieldsManager manager(mesh);

  const int fiberDim = 3;

  const char* labelA = "field A";
  manager.addReal(labelA);
  manager.setFiberDimension(labelA, fiberDim, "vertices");
  const ALE::Obj<real_section_type>& fieldA = manager.getReal(labelA);
  CPPUNIT_ASSERT(!fieldA.isNull());
  const ALE::Obj<Mesh::label_sequence>& vertices = mesh->depthStratum(0);
  mesh->allocate(fieldA);
  for (Mesh::label_sequence::iterator v_iter = vertices->begin();
       v_iter != vertices->end();
       ++v_iter)
    CPPUNIT_ASSERT_EQUAL(fiberDim, fieldA->getFiberDimension(*v_iter));

  const char* labelB = "field B";
  manager.addReal(labelB);
  manager.setFiberDimension(labelB, fiberDim, "cells");
  const ALE::Obj<real_section_type>& fieldB = manager.getReal(labelB);
  CPPUNIT_ASSERT(!fieldB.isNull());
  const ALE::Obj<Mesh::label_sequence>& cells = mesh->heightStratum(0);
  mesh->allocate(fieldB);
  for (Mesh::label_sequence::iterator c_iter = cells->begin();
       c_iter != cells->end();
       ++c_iter)
    CPPUNIT_ASSERT_EQUAL(fiberDim, fieldB->getFiberDimension(*c_iter));
} // testSetFiberDimension

// ----------------------------------------------------------------------
// Test allocate().
void
pylith::topology::TestFieldsManager::testAllocate(void)
{ // testAllocate
  ALE::Obj<Mesh> mesh;
  _initialize(&mesh);
  FieldsManager manager(mesh);

  const int fiberDim = 3;

  const char* labelA = "field A";
  manager.addReal(labelA);
  const ALE::Obj<real_section_type>& fieldA = manager.getReal(labelA);
  manager.setFiberDimension(labelA, fiberDim, "vertices");
  manager.allocate(labelA);

  const char* labelB = "field B";
  manager.addReal(labelB);
  const ALE::Obj<real_section_type>& fieldB = manager.getReal(labelB);
  manager.setFiberDimension(labelB, fiberDim, "cells");
  manager.allocate(labelB);
} // testAllocate

// ----------------------------------------------------------------------
// Test copyLayout().
void
pylith::topology::TestFieldsManager::testCopyLayout(void)
{ // testCopyLayout
  ALE::Obj<Mesh> mesh;
  _initialize(&mesh);
  FieldsManager manager(mesh);
  const int numCells = mesh->heightStratum(0)->size();
  const int offset = numCells;

  const char* labelA = "field A";
  manager.addReal(labelA);
  const int fiberDim = 3;
  const int fixedDim = 1;
  const Mesh::point_type fixedPt = offset + 2;

  const ALE::Obj<real_section_type>& fieldA = manager.getReal(labelA);
  const ALE::Obj<Mesh::label_sequence>& vertices = mesh->depthStratum(0);
  fieldA->setFiberDimension(vertices, fiberDim);
  fieldA->setConstraintDimension(fixedPt, 1);
  mesh->allocate(fieldA);
  fieldA->setConstraintDof(fixedPt, &fixedDim);

  const char* labelB = "field B";
  manager.addReal(labelB);

  manager.copyLayout(labelA);

  size_t size = 2;
  CPPUNIT_ASSERT_EQUAL(size, manager._real.size());
  const ALE::Obj<real_section_type>& fieldB = manager.getReal(labelB);
  
  CPPUNIT_ASSERT_EQUAL(fieldA->size(), fieldB->size());
  CPPUNIT_ASSERT_EQUAL(fieldA->sizeWithBC(), fieldB->sizeWithBC());
  for (Mesh::label_sequence::iterator v_iter = vertices->begin();
       v_iter != vertices->end();
       ++v_iter) {
    if (*v_iter != fixedPt) {
      CPPUNIT_ASSERT_EQUAL(fiberDim, fieldB->getFiberDimension(*v_iter));
      CPPUNIT_ASSERT_EQUAL(0, fieldB->getConstraintDimension(*v_iter));
    } else {
      CPPUNIT_ASSERT_EQUAL(fiberDim, fieldB->getFiberDimension(*v_iter));
      CPPUNIT_ASSERT_EQUAL(1, fieldB->getConstraintDimension(*v_iter));
      CPPUNIT_ASSERT_EQUAL(fixedDim, fieldB->getConstraintDof(*v_iter)[0]);
    } // if/else
  } // for
} // testCopyLayout

// ----------------------------------------------------------------------
// Test copyLayoutFromField().
void
pylith::topology::TestFieldsManager::testCopyLayoutFromField(void)
{ // testCopyLayoutFromField
  ALE::Obj<Mesh> mesh;
  _initialize(&mesh);
  FieldsManager manager(mesh);
  const int numCells = mesh->heightStratum(0)->size();
  const int offset = numCells;

  const char* labelA = "field A";
  const ALE::Obj<real_section_type>& fieldA = 
    new real_section_type(mesh->comm(), mesh->debug());
  const int fiberDim = 3;
  const int fixedDim = 1;
  const Mesh::point_type fixedPt = offset + 2;

  const ALE::Obj<Mesh::label_sequence>& vertices = mesh->depthStratum(0);
  fieldA->setFiberDimension(vertices, fiberDim);
  fieldA->setConstraintDimension(fixedPt, 1);
  mesh->allocate(fieldA);
  fieldA->setConstraintDof(fixedPt, &fixedDim);

  const char* labelB = "field B";
  manager.addReal(labelB);
  const char* labelC = "field C";
  manager.addReal(labelC);

  manager.copyLayout(fieldA);

  size_t size = 2;
  CPPUNIT_ASSERT_EQUAL(size, manager._real.size());
  const ALE::Obj<real_section_type>& fieldB = manager.getReal(labelB);
  const ALE::Obj<real_section_type>& fieldC = manager.getReal(labelC);
  
  CPPUNIT_ASSERT_EQUAL(fieldA->size(), fieldB->size());
  CPPUNIT_ASSERT_EQUAL(fieldA->sizeWithBC(), fieldB->sizeWithBC());
  for (Mesh::label_sequence::iterator v_iter = vertices->begin();
       v_iter != vertices->end();
       ++v_iter) {
    if (*v_iter != fixedPt) {
      CPPUNIT_ASSERT_EQUAL(fiberDim, fieldB->getFiberDimension(*v_iter));
      CPPUNIT_ASSERT_EQUAL(0, fieldB->getConstraintDimension(*v_iter));
      CPPUNIT_ASSERT_EQUAL(fiberDim, fieldC->getFiberDimension(*v_iter));
      CPPUNIT_ASSERT_EQUAL(0, fieldC->getConstraintDimension(*v_iter));
    } else {
      CPPUNIT_ASSERT_EQUAL(fiberDim, fieldB->getFiberDimension(*v_iter));
      CPPUNIT_ASSERT_EQUAL(1, fieldB->getConstraintDimension(*v_iter));
      CPPUNIT_ASSERT_EQUAL(fixedDim, fieldB->getConstraintDof(*v_iter)[0]);
      CPPUNIT_ASSERT_EQUAL(fiberDim, fieldC->getFiberDimension(*v_iter));
      CPPUNIT_ASSERT_EQUAL(1, fieldC->getConstraintDimension(*v_iter));
      CPPUNIT_ASSERT_EQUAL(fixedDim, fieldC->getConstraintDof(*v_iter)[0]);
    } // if/else
  } // for
} // testCopyLayoutFromField

// ----------------------------------------------------------------------
// Test outputField().
void
pylith::topology::TestFieldsManager::testOutputField(void)
{ // testOutputField
  ALE::Obj<Mesh> mesh;
  _initialize(&mesh);
  FieldsManager manager(mesh);

  const std::string& name = "my solution";
  manager.addReal(name.c_str());
  manager.outputField(name.c_str());
  CPPUNIT_ASSERT_EQUAL(name, manager._outputName);
} // testSolutionField

// ----------------------------------------------------------------------
// Test getOutputSoln().
void
pylith::topology::TestFieldsManager::testGetOutputSoln(void)
{ // testGetOutputSoln
  ALE::Obj<Mesh> mesh;
  _initialize(&mesh);
  FieldsManager manager(mesh);

  const char* labels[] = { "field A", "field B", "field C" };
  const int size = 3;
  const int fiberDimA = 2;
  const int fiberDimB = 3;
  const int fiberDimC = 4;

  for (int i=0; i < size; ++i)
    manager.addReal(labels[i]);

  const ALE::Obj<Mesh::label_sequence>& vertices = mesh->depthStratum(0);
  const ALE::Obj<real_section_type>& fieldA = manager.getReal(labels[0]);
  const ALE::Obj<real_section_type>& fieldB = manager.getReal(labels[1]);
  const ALE::Obj<real_section_type>& fieldC = manager.getReal(labels[2]);
  fieldA->setFiberDimension(vertices, fiberDimA);
  fieldB->setFiberDimension(vertices, fiberDimB);
  fieldC->setFiberDimension(vertices, fiberDimC);

  manager.outputField(labels[1]);
  const ALE::Obj<real_section_type>& solution = manager.getOutputSoln();
  CPPUNIT_ASSERT_EQUAL(fiberDimB, 
		       solution->getFiberDimension(*(vertices->begin())));
} // testGetOutputSoln

// ----------------------------------------------------------------------
// Test createHistory().
void
pylith::topology::TestFieldsManager::testCreateHistory(void)
{ // testCreateHistory
  ALE::Obj<Mesh> mesh;
  _initialize(&mesh);
  FieldsManager manager(mesh);

  const char* labels[] = { "field A", "field B", "field C" };
  const int totalSize = 3;
  const int historySize = 2;

  // Add fields
  for (int i=0; i < totalSize; ++i)
    manager.addReal(labels[i]);

  manager.createHistory(labels, historySize);
  for (int i=0; i < historySize; ++i)
    CPPUNIT_ASSERT_EQUAL(std::string(labels[i]), manager._history[i]);
} // testCreateHistory

// ----------------------------------------------------------------------
// Test shiftHistory().
void
pylith::topology::TestFieldsManager::testShiftHistory(void)
{ // testShiftHistory
  ALE::Obj<Mesh> mesh;
  _initialize(&mesh);
  FieldsManager manager(mesh);

  const char* labels[] = { "field A", "field B" };
  const int size = 2;
  const int fiberDimA = 2;
  const int fiberDimB = 3;

  for (int i=0; i < size; ++i)
    manager.addReal(labels[i]);
  manager.createHistory(labels, size);

  const ALE::Obj<Mesh::label_sequence>& vertices = mesh->depthStratum(0);
  const ALE::Obj<real_section_type>& fieldA = manager.getReal(labels[0]);
  const ALE::Obj<real_section_type>& fieldB = manager.getReal(labels[1]);
  fieldA->setFiberDimension(vertices, fiberDimA);
  fieldB->setFiberDimension(vertices, fiberDimB);

  manager.shiftHistory();
  const ALE::Obj<real_section_type>& testA = manager.getReal(labels[0]);
  const ALE::Obj<real_section_type>& testB = manager.getReal(labels[1]);
  CPPUNIT_ASSERT_EQUAL(fiberDimB, 
		       testA->getFiberDimension(*(vertices->begin())));
  CPPUNIT_ASSERT_EQUAL(fiberDimA, 
		       testB->getFiberDimension(*(vertices->begin())));
} // testShiftHistory

// ----------------------------------------------------------------------
// Test getHistoryItem().
void
pylith::topology::TestFieldsManager::testGetHistoryItem(void)
{ // testGetHistoryItem
  ALE::Obj<Mesh> mesh;
  _initialize(&mesh);
  FieldsManager manager(mesh);

  const char* labels[] = { "field A", "field B" };
  const int size = 2;
  const int fiberDimA = 2;
  const int fiberDimB = 3;

  for (int i=0; i < size; ++i)
    manager.addReal(labels[i]);
  manager.createHistory(labels, size);

  const ALE::Obj<Mesh::label_sequence>& vertices = mesh->depthStratum(0);
  const ALE::Obj<real_section_type>& fieldA = manager.getReal(labels[0]);
  const ALE::Obj<real_section_type>& fieldB = manager.getReal(labels[1]);
  fieldA->setFiberDimension(vertices, fiberDimA);
  fieldB->setFiberDimension(vertices, fiberDimB);

  const ALE::Obj<real_section_type>& testA = manager.getHistoryItem(0);
  CPPUNIT_ASSERT_EQUAL(fiberDimA, 
		       testA->getFiberDimension(*(vertices->begin())));

  const ALE::Obj<real_section_type>& testB = manager.getHistoryItem(1);
  CPPUNIT_ASSERT_EQUAL(fiberDimB, 
		       testB->getFiberDimension(*(vertices->begin())));
} // testGetHistoryItem

// ----------------------------------------------------------------------
void
pylith::topology::TestFieldsManager::_initialize(ALE::Obj<Mesh>* mesh) const
{ // _initialize
  CPPUNIT_ASSERT(0 != mesh);

  meshio::MeshIOAscii iohandler;
  iohandler.filename("data/tri3.mesh");
  iohandler.read(mesh);
  CPPUNIT_ASSERT(!mesh->isNull());
} // _initialize


// End of file 
