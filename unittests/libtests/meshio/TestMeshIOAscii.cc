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

#include "TestMeshIOAscii.hh" // Implementation of class methods

#include "pylith/meshio/MeshIOAscii.hh"

// ----------------------------------------------------------------------
CPPUNIT_TEST_SUITE_REGISTRATION( pylith::meshio::TestMeshIOAscii );

// ----------------------------------------------------------------------
// Test constructor
void
pylith::meshio::TestMeshIOAscii::testConstructor(void)
{ // testConstructor
  MeshIOAscii iohandler;
} // testConstructor

// ----------------------------------------------------------------------
// Test filename()
void
pylith::meshio::TestMeshIOAscii::testFilename(void)
{ // testFilename
  MeshIOAscii iohandler;

  const char* filename = "hi.txt";
  iohandler.filename(filename);
  CPPUNIT_ASSERT(0 == strcasecmp(filename, iohandler.filename()));
} // testFilename

// ----------------------------------------------------------------------
// Test write() and read() for 1D mesh.
void
pylith::meshio::TestMeshIOAscii::testWriteRead1D(void)
{ // testWriteRead1D
  CPPUNIT_ASSERT(false);
} // testWriteRead1D

// ----------------------------------------------------------------------
// Test write() and read() for 1D mesh in 2D space.
void
pylith::meshio::TestMeshIOAscii::testWriteRead1Din2D(void)
{ // testWriteRead1D
  CPPUNIT_ASSERT(false);
} // testWriteRead1D

// ----------------------------------------------------------------------
// Test write() and read() for 1D mesh in 3D space.
void
pylith::meshio::TestMeshIOAscii::testWriteRead1Din3D(void)
{ // testWriteRead1D
  CPPUNIT_ASSERT(false);
} // testWriteRead1D

// ----------------------------------------------------------------------
// Test write() and read() for 2D mesh in 2D space.
void
pylith::meshio::TestMeshIOAscii::testWriteRead2D(void)
{ // testWriteRead2D
  CPPUNIT_ASSERT(false);
} // testWriteRead2D

// ----------------------------------------------------------------------
// Test write() and read() for 2D mesh in 3D space.
void
pylith::meshio::TestMeshIOAscii::testWriteRead2Din3D(void)
{ // testWriteRead2Din3D
  CPPUNIT_ASSERT(false);
} // testWriteRead2Din3D

// ----------------------------------------------------------------------
// Test write() and read() for 3D mesh.
void
pylith::meshio::TestMeshIOAscii::testWriteRead3D(void)
{ // testWriteRead3D
  CPPUNIT_ASSERT(false);
} // testWriteRead3D

// End of file 
