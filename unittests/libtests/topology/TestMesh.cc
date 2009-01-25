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

#include "TestMesh.hh" // Implementation of class methods
#include "pylith/topology/Mesh.hh" // USES Mesh

#include "spatialdata/geocoords/CSCart.hh" // USES CSCart

// ----------------------------------------------------------------------
CPPUNIT_TEST_SUITE_REGISTRATION( pylith::topology::TestMesh );

// ----------------------------------------------------------------------
// Test constructor.
void
pylith::topology::TestMesh::testConstructor(void)
{ // testConstructor
  Mesh mesh;
  CPPUNIT_ASSERT(mesh._mesh.isNull());
  CPPUNIT_ASSERT_EQUAL(0, mesh.dimension());
  CPPUNIT_ASSERT_EQUAL(false, mesh.debug());
  CPPUNIT_ASSERT_EQUAL(PETSC_COMM_WORLD, mesh.comm());
  
  Mesh mesh2(2);
  CPPUNIT_ASSERT(!mesh2._mesh.isNull());
  CPPUNIT_ASSERT_EQUAL(2, mesh2.dimension());
  CPPUNIT_ASSERT_EQUAL(PETSC_COMM_WORLD, mesh2.comm());

  Mesh mesh3(1, PETSC_COMM_SELF);
  CPPUNIT_ASSERT(!mesh3._mesh.isNull());
  CPPUNIT_ASSERT_EQUAL(1, mesh3.dimension());
  CPPUNIT_ASSERT_EQUAL(PETSC_COMM_SELF, mesh3.comm());
} // testConstructor

// ----------------------------------------------------------------------
// Test constructor.
void
pylith::topology::TestMesh::testCreateSieveMesh(void)
{ // testCreateSieveMesh
  Mesh mesh;
  CPPUNIT_ASSERT(mesh._mesh.isNull());

  int dim = 2;
  mesh.createSieveMesh(dim);
  CPPUNIT_ASSERT(!mesh._mesh.isNull());
  CPPUNIT_ASSERT_EQUAL(dim, mesh.dimension());

  dim = 1;
  mesh.createSieveMesh(dim);
  CPPUNIT_ASSERT(!mesh._mesh.isNull());
  CPPUNIT_ASSERT_EQUAL(dim, mesh.dimension());
} // testCreateMeshSieve

// ----------------------------------------------------------------------
// Test sieveMesh().
void
pylith::topology::TestMesh::testSieveMesh(void)
{ // testSieveMesh
  const int dim = 2;
  Mesh mesh(dim);
  
  const ALE::Obj<SieveMesh>& sieveMesh = mesh.sieveMesh();
  CPPUNIT_ASSERT(!sieveMesh.isNull());
  CPPUNIT_ASSERT_EQUAL(dim, sieveMesh->getDimension());
} // testSieveMesh

// ----------------------------------------------------------------------
// Test coordsys().
void
pylith::topology::TestMesh::testCoordsys(void)
{ // testCoordsys
  Mesh mesh;

  spatialdata::geocoords::CSCart cs;
  cs.setSpaceDim(2);

  mesh.coordsys(&cs);

  CPPUNIT_ASSERT_EQUAL(cs.spaceDim(), mesh.coordsys()->spaceDim());
} // testCoordsys

// ----------------------------------------------------------------------
// Test debug().
void
pylith::topology::TestMesh::testDebug(void)
{ // testDebug
  Mesh mesh;
  CPPUNIT_ASSERT_EQUAL(false, mesh.debug());

  mesh.debug(true);
  CPPUNIT_ASSERT_EQUAL(true, mesh.debug());
} // testDebug

// ----------------------------------------------------------------------
// Test dimension().
void
pylith::topology::TestMesh::testDimension(void)
{ // testDimension
  Mesh mesh;
  CPPUNIT_ASSERT_EQUAL(0, mesh.dimension());

  const int dim = 2;
  Mesh mesh2(dim);
  CPPUNIT_ASSERT_EQUAL(dim, mesh2.dimension());
} // testDimension

// ----------------------------------------------------------------------
// Test comm().
void
pylith::topology::TestMesh::testComm(void)
{ // testComm
  Mesh mesh;
  CPPUNIT_ASSERT_EQUAL(PETSC_COMM_WORLD, mesh.comm());

  mesh.comm(PETSC_COMM_SELF);
  CPPUNIT_ASSERT_EQUAL(PETSC_COMM_SELF, mesh.comm());
} // testComm

// ----------------------------------------------------------------------
// Test initialize().
void
pylith::topology::TestMesh::testInitialize(void)
{ // testInitialize
  Mesh mesh;

  spatialdata::geocoords::CSCart cs;
  cs.setSpaceDim(2);

  mesh.coordsys(&cs);
  mesh.initialize();

} // testInitialize


// End of file 
