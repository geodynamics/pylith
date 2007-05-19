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

#include <stdexcept> // TEMPORARY

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
// Test copyLayout().
void
pylith::topology::TestFieldsManager::testCopyLayout(void)
{ // testCopyLayout
  throw std::logic_error("Unit test not implemented.");
} // testCopyLayout

// ----------------------------------------------------------------------
// Test copyLayoutFromField().
void
pylith::topology::TestFieldsManager::testCopyLayoutFromField(void)
{ // testCopyLayoutFromField
  throw std::logic_error("Unit test not implemented.");
} // testCopyLayoutFromField

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
