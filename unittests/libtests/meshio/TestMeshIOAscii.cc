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

#include "data/MeshData1D.hh"
#include "data/MeshData1Din2D.hh"
#include "data/MeshData1Din3D.hh"
#include "data/MeshData2D.hh"
#include "data/MeshData2Din3D.hh"
#include "data/MeshData3D.hh"

#include <assert.h> // USES assert()

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
  MeshData1D data;
  const char* filename = "mesh1D.txt";
  _testWriteRead(data, filename);
} // testWriteRead1D

// ----------------------------------------------------------------------
// Test write() and read() for 1D mesh in 2D space.
void
pylith::meshio::TestMeshIOAscii::testWriteRead1Din2D(void)
{ // testWriteRead1Din2D
  MeshData1Din2D data;
  const char* filename = "mesh1Din2D.txt";
  _testWriteRead(data, filename);
} // testWriteRead1Din2D

// ----------------------------------------------------------------------
// Test write() and read() for 1D mesh in 3D space.
void
pylith::meshio::TestMeshIOAscii::testWriteRead1Din3D(void)
{ // testWriteRead1Din3D
  MeshData1Din3D data;
  const char* filename = "mesh1Din3D.txt";
  _testWriteRead(data, filename);
} // testWriteRead1Din3D

// ----------------------------------------------------------------------
// Test write() and read() for 2D mesh in 2D space.
void
pylith::meshio::TestMeshIOAscii::testWriteRead2D(void)
{ // testWriteRead2D
  MeshData2D data;
  const char* filename = "mesh2D.txt";
  _testWriteRead(data, filename);
} // testWriteRead2D

// ----------------------------------------------------------------------
// Test write() and read() for 2D mesh in 3D space.
void
pylith::meshio::TestMeshIOAscii::testWriteRead2Din3D(void)
{ // testWriteRead2Din3D
  MeshData2Din3D data;
  const char* filename = "mesh2Din3D.txt";
  _testWriteRead(data, filename);
} // testWriteRead2Din3D

// ----------------------------------------------------------------------
// Test write() and read() for 3D mesh.
void
pylith::meshio::TestMeshIOAscii::testWriteRead3D(void)
{ // testWriteRead3D
  MeshData3D data;
  const char* filename = "mesh3D.txt";
  _testWriteRead(data, filename);
} // testWriteRead3D

// ----------------------------------------------------------------------
// Build mesh, perform write() and read(), and then check values.
void
pylith::meshio::TestMeshIOAscii::_testWriteRead(const MeshData& data,
						const char* filename)
{ // _testWriteRead
  ALE::Obj<Mesh>* meshOut = createMesh(data);

  // Write mesh
  MeshIOAscii iohandler;
  iohandler.filename(filename);
  iohandler.write(meshOut);
  delete meshOut; meshOut = 0;

  // Read mesh
  ALE::Obj<Mesh> meshIn;
  iohandler.read(&meshIn);

  // Make sure meshIn matches data
  checkVals(meshIn, data);
} // _testWriteRead

// End of file 
