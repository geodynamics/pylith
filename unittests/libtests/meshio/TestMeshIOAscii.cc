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

#include "TestMeshIOAscii.hh" // Implementation of class methods

#include "pylith/meshio/MeshIOAscii.hh"

#include "pylith/topology/Mesh.hh" // USES Mesh

#include "data/MeshData1D.hh"
#include "data/MeshData1Din2D.hh"
#include "data/MeshData1Din3D.hh"
#include "data/MeshData2D.hh"
#include "data/MeshData2Din3D.hh"
#include "data/MeshData3D.hh"
#include "data/MeshData3DIndexOne.hh"

#include <strings.h> // USES strcasecmp()

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
// Test debug()
void
pylith::meshio::TestMeshIOAscii::testDebug(void)
{ // testDebug
  MeshIOAscii iohandler;
  _testDebug(iohandler);
} // testDebug

// ----------------------------------------------------------------------
// Test interpolate()
void
pylith::meshio::TestMeshIOAscii::testInterpolate(void)
{ // testInterpolate
  MeshIOAscii iohandler;
  _testInterpolate(iohandler);
} // testInterpolate

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
// Test read() for 3D mesh with 1 based indexing.
void
pylith::meshio::TestMeshIOAscii::testRead3DIndexOne(void)
{ // testRead3IndexDOne
  MeshData3DIndexOne data;
  const char* filename = "data/mesh3DIndexOne.txt";
  _testRead(data, filename);
} // testRead3DIndexOne

// ----------------------------------------------------------------------
// Test write() and read() for 2D mesh in 2D space with comments.
void
pylith::meshio::TestMeshIOAscii::testReadComments(void)
{ // testWriteReadComments
  MeshData2D data;
  const char* filename = "data/mesh2D_comments.txt";
  _testRead(data, filename);
} // testWriteReadComments

// ----------------------------------------------------------------------
// Build mesh, perform write() and read(), and then check values.
void
pylith::meshio::TestMeshIOAscii::_testWriteRead(const MeshData& data,
						const char* filename)
{ // _testWriteRead
  topology::Mesh* meshOut = _createMesh(data);

  // Write mesh
  MeshIOAscii iohandler;
  iohandler.filename(filename);
  iohandler.write(meshOut);
  delete meshOut; meshOut = 0;

  // Read mesh
  topology::Mesh meshIn;
  iohandler.read(&meshIn);

  // Make sure meshIn matches data
  _checkVals(meshIn, data);
} // _testWriteRead

// ----------------------------------------------------------------------
// Read mesh and then check values.
void
pylith::meshio::TestMeshIOAscii::_testRead(const MeshData& data,
					   const char* filename)
{ // _testWriteRead
  // Read mesh
  topology::Mesh mesh;
  MeshIOAscii iohandler;
  iohandler.filename(filename);
  iohandler.read(&mesh);

  // Make sure mesh matches data
  _checkVals(mesh, data);
} // _testRead


// End of file 
