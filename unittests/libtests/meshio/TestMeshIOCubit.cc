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

#include "TestMeshIOCubit.hh" // Implementation of class methods

#include "pylith/meshio/MeshIOCubit.hh"

#include "pylith/utils/sievetypes.hh" // USES PETSc Mesh
#include "pylith/utils/array.hh" // USES int_array

#include "data/MeshDataCubitTri.hh"
#include "data/MeshDataCubitQuad.hh"
#include "data/MeshDataCubitTet.hh"
#include "data/MeshDataCubitHex.hh"

// ----------------------------------------------------------------------
CPPUNIT_TEST_SUITE_REGISTRATION( pylith::meshio::TestMeshIOCubit );

// ----------------------------------------------------------------------
// Test constructor
void
pylith::meshio::TestMeshIOCubit::testConstructor(void)
{ // testConstructor
  MeshIOCubit iohandler;
} // testConstructor

// ----------------------------------------------------------------------
// Test debug()
void
pylith::meshio::TestMeshIOCubit::testDebug(void)
{ // testDebug
  MeshIOCubit iohandler;
  _testDebug(iohandler);
} // testDebug

// ----------------------------------------------------------------------
// Test interpolate()
void
pylith::meshio::TestMeshIOCubit::testInterpolate(void)
{ // testInterpolate
  MeshIOCubit iohandler;
  _testInterpolate(iohandler);
} // testInterpolate

// ----------------------------------------------------------------------
// Test filename()
void
pylith::meshio::TestMeshIOCubit::testFilename(void)
{ // testFilename
  MeshIOCubit iohandler;

  const char* filename = "hi.txt";
  iohandler.filename(filename);
  CPPUNIT_ASSERT(0 == strcasecmp(filename, iohandler.filename()));
} // testFilename

// ----------------------------------------------------------------------
// Test read() for mesh with triangle cells.
void
pylith::meshio::TestMeshIOCubit::testReadTri(void)
{ // testReadTri
  MeshDataCubitTri data;
  const char* filename = "data/twotri3.exo";
  _testRead(data, filename);
} // testReadTri

// ----------------------------------------------------------------------
// Test read() for mesh with quadrilateral cells.
void
pylith::meshio::TestMeshIOCubit::testReadQuad(void)
{ // testReadQuad
  MeshDataCubitQuad data;
  const char* filename = "data/twoquad4.exo";
  _testRead(data, filename);
} // testReadQuad

// ----------------------------------------------------------------------
// Test read() for mesh with tetrahedral cells.
void
pylith::meshio::TestMeshIOCubit::testReadTet(void)
{ // testReadTet
  MeshDataCubitTet data;
  const char* filename = "data/twotet4.exo";
  _testRead(data, filename);
} // testReadTet

// ----------------------------------------------------------------------
// Test read() for mesh with tetrahedral cells.
void
pylith::meshio::TestMeshIOCubit::testReadHex(void)
{ // testReadHex
  MeshDataCubitHex data;
  const char* filename = "data/twohex8.exo";
  _testRead(data, filename);
} // testReadHex

// ----------------------------------------------------------------------
// Build mesh, perform read(), and then check values.
void
pylith::meshio::TestMeshIOCubit::_testRead(const MeshData& data,
					   const char* filename)
{ // _testRead
  MeshIOCubit iohandler;
  iohandler.filename(filename);

  // Read mesh
  ALE::Obj<Mesh> mesh;
  iohandler.read(&mesh);

  // Make sure meshIn matches data
  _checkVals(mesh, data);
} // _testRead

// ----------------------------------------------------------------------
// Test _orientCells with line cells.
void
pylith::meshio::TestMeshIOCubit::testOrientLine(void)
{ // testOrientLine
  // No change in cells exepected

  const int meshDim = 1;
  const int numCells = 2;
  const int numCorners = 2;
  const int cellsOrig[] = {
    0, 1,
    2, 3
  };
  const int cellsE[] = {
    0, 1,
    2, 3
  };
  
  int_array cells(cellsOrig, numCells*numCorners);
  MeshIOCubit::_orientCells(&cells, numCells, numCorners, meshDim);

  const int size = numCells*numCorners;
  CPPUNIT_ASSERT_EQUAL(size, int(cells.size()));
  for (int i=0; i < size; ++i)
    CPPUNIT_ASSERT_EQUAL(cellsE[i], cells[i]);
} // testOrientLine

// ----------------------------------------------------------------------
// Test _orientCells with tri cells.
void
pylith::meshio::TestMeshIOCubit::testOrientTri(void)
{ // testOrientTri
  // Swap vertices 1 and 2, vertex 0 is unchanged.

  const int meshDim = 2;
  const int numCells = 2;
  const int numCorners = 3;
  const int cellsOrig[] = {
    0, 1, 2,
    3, 4, 5
  };
  const int cellsE[] = {
    0, 2, 1,
    3, 5, 4
  };
  
  int_array cells(cellsOrig, numCells*numCorners);
  MeshIOCubit::_orientCells(&cells, numCells, numCorners, meshDim);

  const int size = numCells*numCorners;
  CPPUNIT_ASSERT_EQUAL(size, int(cells.size()));
  for (int i=0; i < size; ++i)
    CPPUNIT_ASSERT_EQUAL(cellsE[i], cells[i]);
} // testOrientTri

// ----------------------------------------------------------------------
// Test _orientCells with quad cells.
void
pylith::meshio::TestMeshIOCubit::testOrientQuad(void)
{ // testOrientQuad
  // Expect vertices 0 and 1 to remain the same, but vertices 2 and 3
  // should be swapped

  const int meshDim = 2;
  const int numCells = 2;
  const int numCorners = 4;
  const int cellsOrig[] = {
    0, 1, 2, 3,
    6, 7, 8, 9
  };
  const int cellsE[] = {
    0, 1, 3, 2,
    6, 7, 9, 8
  };
  
  int_array cells(cellsOrig, numCells*numCorners);
  MeshIOCubit::_orientCells(&cells, numCells, numCorners, meshDim);

  const int size = numCells*numCorners;
  CPPUNIT_ASSERT_EQUAL(size, int(cells.size()));
  for (int i=0; i < size; ++i)
    CPPUNIT_ASSERT_EQUAL(cellsE[i], cells[i]);
} // testOrientQuad

// ----------------------------------------------------------------------
// Test _orientCells with tet cells.
void
pylith::meshio::TestMeshIOCubit::testOrientTet(void)
{ // testOrientTet
  // No change in cells exepected

  const int meshDim = 3;
  const int numCells = 2;
  const int numCorners = 4;
  const int cellsOrig[] = {
    0, 1, 2, 3,
    3, 4, 5, 6
  };
  const int cellsE[] = {
    0, 1, 2, 3,
    3, 4, 5, 6
  };
  
  int_array cells(cellsOrig, numCells*numCorners);
  MeshIOCubit::_orientCells(&cells, numCells, numCorners, meshDim);

  const int size = numCells*numCorners;
  CPPUNIT_ASSERT_EQUAL(size, int(cells.size()));
  for (int i=0; i < size; ++i)
    CPPUNIT_ASSERT_EQUAL(cellsE[i], cells[i]);
} // testOrientTet

// ----------------------------------------------------------------------
// Test _orientCells with hex cells.
void
pylith::meshio::TestMeshIOCubit::testOrientHex(void)
{ // testOrientHex
  // Expect vertices 0 and 1 to remain the same, but vertices 2 and 3
  // should be swapped

  const int meshDim = 3;
  const int numCells = 2;
  const int numCorners = 8;
  const int cellsOrig[] = {
    0, 1, 2, 3, 4, 5, 6, 7,
    10, 11, 12, 13, 14, 15, 16, 17
  };
  const int cellsE[] = {
    0, 1, 3, 2, 4, 5, 7, 6,
    10, 11, 13, 12, 14, 15, 17, 16
  };
  
  int_array cells(cellsOrig, numCells*numCorners);
  MeshIOCubit::_orientCells(&cells, numCells, numCorners, meshDim);

  const int size = numCells*numCorners;
  CPPUNIT_ASSERT_EQUAL(size, int(cells.size()));
  for (int i=0; i < size; ++i)
    CPPUNIT_ASSERT_EQUAL(cellsE[i], cells[i]);
} // testOrientHex


// End of file 
