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

#include "TestMeshIOCubit.hh" // Implementation of class methods

#include "pylith/meshio/MeshIOCubit.hh"

#include "pylith/topology/Mesh.hh" // USES Mesh
#include "pylith/utils/array.hh" // USES int_array

#include "data/MeshDataCubitTri.hh"
#include "data/MeshDataCubitQuad.hh"
#include "data/MeshDataCubitTet.hh"
#include "data/MeshDataCubitHex.hh"

#include <strings.h> // USES strcasecmp()

// ----------------------------------------------------------------------
CPPUNIT_TEST_SUITE_REGISTRATION( pylith::meshio::TestMeshIOCubit );

// ----------------------------------------------------------------------
// Test constructor
void
pylith::meshio::TestMeshIOCubit::testConstructor(void)
{ // testConstructor
  PYLITH_METHOD_BEGIN;

  MeshIOCubit iohandler;

  PYLITH_METHOD_END;
} // testConstructor

// ----------------------------------------------------------------------
// Test debug()
void
pylith::meshio::TestMeshIOCubit::testDebug(void)
{ // testDebug
  PYLITH_METHOD_BEGIN;

  MeshIOCubit iohandler;
  _testDebug(iohandler);

  PYLITH_METHOD_END;
} // testDebug

// ----------------------------------------------------------------------
// Test interpolate()
void
pylith::meshio::TestMeshIOCubit::testInterpolate(void)
{ // testInterpolate
  PYLITH_METHOD_BEGIN;

  MeshIOCubit iohandler;
  _testInterpolate(iohandler);

  PYLITH_METHOD_END;
} // testInterpolate

// ----------------------------------------------------------------------
// Test filename()
void
pylith::meshio::TestMeshIOCubit::testFilename(void)
{ // testFilename
  PYLITH_METHOD_BEGIN;

  MeshIOCubit iohandler;

  const char* filename = "hi.txt";
  iohandler.filename(filename);
  CPPUNIT_ASSERT(0 == strcasecmp(filename, iohandler.filename()));

  PYLITH_METHOD_END;
} // testFilename

// ----------------------------------------------------------------------
// Test read() for mesh with triangle cells.
void
pylith::meshio::TestMeshIOCubit::testReadTri(void)
{ // testReadTri
  PYLITH_METHOD_BEGIN;

  MeshDataCubitTri data;
  _testRead(data, "data/twotri3_12.2.exo");
  _testRead(data, "data/twotri3_13.0.exo");

  PYLITH_METHOD_END;
} // testReadTri

// ----------------------------------------------------------------------
// Test read() for mesh with quadrilateral cells.
void
pylith::meshio::TestMeshIOCubit::testReadQuad(void)
{ // testReadQuad
  PYLITH_METHOD_BEGIN;

  MeshDataCubitQuad data;
  _testRead(data, "data/twoquad4_12.2.exo");
  _testRead(data, "data/twoquad4_13.0.exo");

  PYLITH_METHOD_END;
} // testReadQuad

// ----------------------------------------------------------------------
// Test read() for mesh with tetrahedral cells.
void
pylith::meshio::TestMeshIOCubit::testReadTet(void)
{ // testReadTet
  PYLITH_METHOD_BEGIN;

  MeshDataCubitTet data;
  _testRead(data, "data/twotet4_12.2.exo");
  _testRead(data, "data/twotet4_13.0.exo");

  PYLITH_METHOD_END;
} // testReadTet

// ----------------------------------------------------------------------
// Test read() for mesh with tetrahedral cells.
void
pylith::meshio::TestMeshIOCubit::testReadHex(void)
{ // testReadHex
  PYLITH_METHOD_BEGIN;

  MeshDataCubitHex data;
  _testRead(data, "data/twohex8_12.2.exo");
  _testRead(data, "data/twohex8_13.0.exo");

  PYLITH_METHOD_END;
} // testReadHex

// ----------------------------------------------------------------------
// Build mesh, perform read(), and then check values.
void
pylith::meshio::TestMeshIOCubit::_testRead(const MeshData& data,
					   const char* filename)
{ // _testRead
  PYLITH_METHOD_BEGIN;

  MeshIOCubit iohandler;
  iohandler.filename(filename);
  iohandler.useNodesetNames(true);

  // Read mesh
  delete _mesh; _mesh = new topology::Mesh;
  iohandler.read(_mesh);

  // Make sure mesh matches data
  _checkVals(data);

  PYLITH_METHOD_END;
} // _testRead

// ----------------------------------------------------------------------
// Test _orientCells with line cells.
void
pylith::meshio::TestMeshIOCubit::testOrientLine(void)
{ // testOrientLine
  PYLITH_METHOD_BEGIN;

  // No change in cells exepected

  const int meshDim = 1;
  const int numCells = 2;
  const int numCorners = 2;
  const int cellsOrig[numCells*numCorners] = {
    0, 1,
    2, 3
  };
  const int cellsE[numCells*numCorners] = {
    0, 1,
    2, 3
  };
  
  int_array cells(cellsOrig, numCells*numCorners);
  MeshIOCubit::_orientCells(&cells, numCells, numCorners, meshDim);

  const int size = numCells*numCorners;
  CPPUNIT_ASSERT_EQUAL(size, int(cells.size()));
  for (int i=0; i < size; ++i)
    CPPUNIT_ASSERT_EQUAL(cellsE[i], cells[i]);

  PYLITH_METHOD_END;
} // testOrientLine

// ----------------------------------------------------------------------
// Test _orientCells with tri cells.
void
pylith::meshio::TestMeshIOCubit::testOrientTri(void)
{ // testOrientTri
  PYLITH_METHOD_BEGIN;

  // No changes

  const int meshDim = 2;
  const int numCells = 2;
  const int numCorners = 3;
  const int cellsOrig[numCells*numCorners] = {
    0, 1, 2,
    3, 4, 5
  };
  const int cellsE[numCells*numCorners] = {
    0, 1, 2,
    3, 4, 5
  };
  
  int_array cells(cellsOrig, numCells*numCorners);
  MeshIOCubit::_orientCells(&cells, numCells, numCorners, meshDim);

  const int size = numCells*numCorners;
  CPPUNIT_ASSERT_EQUAL(size, int(cells.size()));
  for (int i=0; i < size; ++i)
    CPPUNIT_ASSERT_EQUAL(cellsE[i], cells[i]);

  PYLITH_METHOD_END;
} // testOrientTri

// ----------------------------------------------------------------------
// Test _orientCells with quad cells.
void
pylith::meshio::TestMeshIOCubit::testOrientQuad(void)
{ // testOrientQuad
  PYLITH_METHOD_BEGIN;

  // Expect no change.

  const int meshDim = 2;
  const int numCells = 2;
  const int numCorners = 4;
  const int cellsOrig[numCells*numCorners] = {
    0, 1, 2, 3,
    6, 7, 8, 9
  };
  const int cellsE[numCells*numCorners] = {
    0, 1, 2, 3,
    6, 7, 8, 9
  };
  
  int_array cells(cellsOrig, numCells*numCorners);
  MeshIOCubit::_orientCells(&cells, numCells, numCorners, meshDim);

  const int size = numCells*numCorners;
  CPPUNIT_ASSERT_EQUAL(size, int(cells.size()));
  for (int i=0; i < size; ++i)
    CPPUNIT_ASSERT_EQUAL(cellsE[i], cells[i]);

  PYLITH_METHOD_END;
} // testOrientQuad

// ----------------------------------------------------------------------
// Test _orientCells with tet cells.
void
pylith::meshio::TestMeshIOCubit::testOrientTet(void)
{ // testOrientTet
  PYLITH_METHOD_BEGIN;

  // No change in cells exepected

  const int meshDim = 3;
  const int numCells = 2;
  const int numCorners = 4;
  const int cellsOrig[numCells*numCorners] = {
    0, 1, 2, 3,
    3, 4, 5, 6
  };
  const int cellsE[numCells*numCorners] = {
    0, 1, 2, 3,
    3, 4, 5, 6
  };
  
  int_array cells(cellsOrig, numCells*numCorners);
  MeshIOCubit::_orientCells(&cells, numCells, numCorners, meshDim);

  const int size = numCells*numCorners;
  CPPUNIT_ASSERT_EQUAL(size, int(cells.size()));
  for (int i=0; i < size; ++i)
    CPPUNIT_ASSERT_EQUAL(cellsE[i], cells[i]);

  PYLITH_METHOD_END;
} // testOrientTet

// ----------------------------------------------------------------------
// Test _orientCells with hex cells.
void
pylith::meshio::TestMeshIOCubit::testOrientHex(void)
{ // testOrientHex
  PYLITH_METHOD_BEGIN;

  // Expect no change.

  const int meshDim = 3;
  const int numCells = 2;
  const int numCorners = 8;
  const int cellsOrig[] = {
    0, 1, 2, 3, 4, 5, 6, 7,
    10, 11, 12, 13, 14, 15, 16, 17
  };
  const int cellsE[] = {
    0, 1, 2, 3, 4, 5, 6, 7,
    10, 11, 12, 13, 14, 15, 16, 17
  };
  
  int_array cells(cellsOrig, numCells*numCorners);
  MeshIOCubit::_orientCells(&cells, numCells, numCorners, meshDim);

  const int size = numCells*numCorners;
  CPPUNIT_ASSERT_EQUAL(size, int(cells.size()));
  for (int i=0; i < size; ++i)
    CPPUNIT_ASSERT_EQUAL(cellsE[i], cells[i]);

  PYLITH_METHOD_END;
} // testOrientHex


// End of file 
