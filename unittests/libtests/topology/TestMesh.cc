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
} // testConstructor

// ----------------------------------------------------------------------
// Test sieveMesh().
void
pylith::topology::TestMesh::testSieveMesh(void)
{ // testSieveMesh
  const int dim = 2;

  Mesh mesh(PETSC_COMM_WORLD, dim);
  
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
