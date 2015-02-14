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
// Copyright (c) 2010-2015 University of California, Davis
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
  PYLITH_METHOD_BEGIN;

  MeshIOAscii iohandler;

  PYLITH_METHOD_END;
} // testConstructor

// ----------------------------------------------------------------------
// Test debug()
void
pylith::meshio::TestMeshIOAscii::testDebug(void)
{ // testDebug
  PYLITH_METHOD_BEGIN;

  MeshIOAscii iohandler;
  _testDebug(iohandler);

  PYLITH_METHOD_END;
} // testDebug

// ----------------------------------------------------------------------
// Test interpolate()
void
pylith::meshio::TestMeshIOAscii::testInterpolate(void)
{ // testInterpolate
  PYLITH_METHOD_BEGIN;

  MeshIOAscii iohandler;
  _testInterpolate(iohandler);

  PYLITH_METHOD_END;
} // testInterpolate

// ----------------------------------------------------------------------
// Test filename()
void
pylith::meshio::TestMeshIOAscii::testFilename(void)
{ // testFilename
  PYLITH_METHOD_BEGIN;

  MeshIOAscii iohandler;

  const char* filename = "hi.txt";
  iohandler.filename(filename);
  CPPUNIT_ASSERT(0 == strcasecmp(filename, iohandler.filename()));

  PYLITH_METHOD_END;
} // testFilename

// ----------------------------------------------------------------------
// Test write() and read() for 1D mesh.
void
pylith::meshio::TestMeshIOAscii::testWriteRead1D(void)
{ // testWriteRead1D
  PYLITH_METHOD_BEGIN;

  MeshData1D data;
  const char* filename = "mesh1D.txt";
  _testWriteRead(data, filename);

  PYLITH_METHOD_END;
} // testWriteRead1D

// ----------------------------------------------------------------------
// Test write() and read() for 1D mesh in 2D space.
void
pylith::meshio::TestMeshIOAscii::testWriteRead1Din2D(void)
{ // testWriteRead1Din2D
  PYLITH_METHOD_BEGIN;

  MeshData1Din2D data;
  const char* filename = "mesh1Din2D.txt";
  _testWriteRead(data, filename);

  PYLITH_METHOD_END;
} // testWriteRead1Din2D

// ----------------------------------------------------------------------
// Test write() and read() for 1D mesh in 3D space.
void
pylith::meshio::TestMeshIOAscii::testWriteRead1Din3D(void)
{ // testWriteRead1Din3D
  PYLITH_METHOD_BEGIN;

  MeshData1Din3D data;
  const char* filename = "mesh1Din3D.txt";
  _testWriteRead(data, filename);

  PYLITH_METHOD_END;
} // testWriteRead1Din3D

// ----------------------------------------------------------------------
// Test write() and read() for 2D mesh in 2D space.
void
pylith::meshio::TestMeshIOAscii::testWriteRead2D(void)
{ // testWriteRead2D
  PYLITH_METHOD_BEGIN;

  MeshData2D data;
  const char* filename = "mesh2D.txt";
  _testWriteRead(data, filename);

  PYLITH_METHOD_END;
} // testWriteRead2D

// ----------------------------------------------------------------------
// Test write() and read() for 2D mesh in 3D space.
void
pylith::meshio::TestMeshIOAscii::testWriteRead2Din3D(void)
{ // testWriteRead2Din3D
  PYLITH_METHOD_BEGIN;

  MeshData2Din3D data;
  const char* filename = "mesh2Din3D.txt";
  _testWriteRead(data, filename);

  PYLITH_METHOD_END;
} // testWriteRead2Din3D

// ----------------------------------------------------------------------
// Test write() and read() for 3D mesh.
void
pylith::meshio::TestMeshIOAscii::testWriteRead3D(void)
{ // testWriteRead3D
  PYLITH_METHOD_BEGIN;

  MeshData3D data;
  const char* filename = "mesh3D.txt";
  _testWriteRead(data, filename);

  PYLITH_METHOD_END;
} // testWriteRead3D

// ----------------------------------------------------------------------
// Test read() for 3D mesh with 1 based indexing.
void
pylith::meshio::TestMeshIOAscii::testRead3DIndexOne(void)
{ // testRead3IndexDOne
  PYLITH_METHOD_BEGIN;

  MeshData3DIndexOne data;
  const char* filename = "data/mesh3DIndexOne.txt";
  _testRead(data, filename);

  PYLITH_METHOD_END;
} // testRead3DIndexOne

// ----------------------------------------------------------------------
// Test write() and read() for 2D mesh in 2D space with comments.
void
pylith::meshio::TestMeshIOAscii::testReadComments(void)
{ // testWriteReadComments
  PYLITH_METHOD_BEGIN;

  MeshData2D data;
  const char* filename = "data/mesh2D_comments.txt";
  _testRead(data, filename);

  PYLITH_METHOD_END;
} // testWriteReadComments

// ----------------------------------------------------------------------
// Build mesh, perform write() and read(), and then check values.
void
pylith::meshio::TestMeshIOAscii::_testWriteRead(const MeshData& data,
						const char* filename)
{ // _testWriteRead
  PYLITH_METHOD_BEGIN;

  _createMesh(data);

  // Write mesh
  MeshIOAscii iohandler;
  iohandler.filename(filename);
  iohandler.write(_mesh);

  // Read mesh
  delete _mesh; _mesh = new topology::Mesh;
  iohandler.read(_mesh);

  // Make sure meshIn matches data
  _checkVals(data);

  PYLITH_METHOD_END;
} // _testWriteRead

// ----------------------------------------------------------------------
// Read mesh and then check values.
void
pylith::meshio::TestMeshIOAscii::_testRead(const MeshData& data,
					   const char* filename)
{ // _testWriteRead
  PYLITH_METHOD_BEGIN;

  // Read mesh
  MeshIOAscii iohandler;
  iohandler.filename(filename);
  delete _mesh; _mesh = new topology::Mesh;
  iohandler.read(_mesh);

  // Make sure mesh matches data
  _checkVals(data);

  PYLITH_METHOD_END;
} // _testRead


// End of file 
