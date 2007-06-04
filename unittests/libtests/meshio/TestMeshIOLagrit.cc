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

#include "TestMeshIOLagrit.hh" // Implementation of class methods

#include "pylith/meshio/MeshIOLagrit.hh"

#include "pylith/utils/sievetypes.hh" // USES PETSc Mesh
#include "pylith/utils/array.hh" // USES int_array

#include "data/MeshDataLagritTet.hh"

#if defined(WORDS_BIGENDIAN)
#define NATIVE_BIG_ENDIAN
#else
#define NATIVE_LITTLE_ENDIAN
#endif

// ----------------------------------------------------------------------
CPPUNIT_TEST_SUITE_REGISTRATION( pylith::meshio::TestMeshIOLagrit );

// ----------------------------------------------------------------------
// Test constructor
void
pylith::meshio::TestMeshIOLagrit::testConstructor(void)
{ // testConstructor
  MeshIOLagrit iohandler;
} // testConstructor

// ----------------------------------------------------------------------
// Test debug()
void
pylith::meshio::TestMeshIOLagrit::testDebug(void)
{ // testDebug
  MeshIOLagrit iohandler;
  _testDebug(iohandler);
} // testDebug

// ----------------------------------------------------------------------
// Test interpolate()
void
pylith::meshio::TestMeshIOLagrit::testInterpolate(void)
{ // testInterpolate
  MeshIOLagrit iohandler;
  _testInterpolate(iohandler);
} // testInterpolate

// ----------------------------------------------------------------------
// Test filename()
void
pylith::meshio::TestMeshIOLagrit::testFilename(void)
{ // testFilename
  MeshIOLagrit iohandler;

  const char* filenameGmv = "hi.txt";
  const char* filenamePset = "hi2.txt";
  iohandler.filenameGmv(filenameGmv);
  iohandler.filenamePset(filenamePset);

  CPPUNIT_ASSERT(0 == strcasecmp(filenameGmv, iohandler.filenameGmv()));
  CPPUNIT_ASSERT(0 == strcasecmp(filenamePset, iohandler.filenamePset()));
} // testFilename

// ----------------------------------------------------------------------
// Test read() for mesh with ASCII files.
void
pylith::meshio::TestMeshIOLagrit::testReadTetAscii(void)
{ // testReadTetAscii
  MeshDataLagritTet data;
  const char* filenameGmv = "data/cube2_ascii.gmv";
  const char* filenamePset = "data/cube2_ascii.pset";
  _testRead(data, filenameGmv, filenamePset);
} // testReadTetAscii

// ----------------------------------------------------------------------
// Test read() for mesh with binary files.
void
pylith::meshio::TestMeshIOLagrit::testReadTetBinary(void)
{ // testReadTetBinary
  MeshDataLagritTet data;
  const char* filenameGmv = "data/cube2_binary.gmv";
  const char* filenamePset = "data/cube2_binary.pset";
  _testRead(data, filenameGmv, filenamePset);
} // testReadTetBinary

// ----------------------------------------------------------------------
// Build mesh, perform read(), and then check values.
void
pylith::meshio::TestMeshIOLagrit::_testRead(const MeshData& data,
					    const char* filenameGmv,
					    const char* filenamePset)
{ // _testRead
  MeshIOLagrit iohandler;
  iohandler.filenameGmv(filenameGmv);
  iohandler.filenamePset(filenamePset);
  
  // LaGriT file was created on little endian machine, so flip endian if
  // running test on big endian machine.
#if defined(NATIVE_LITTLE_ENDIAN)
  iohandler.flipEndian(false);
#else
  iohandler.flipEndian(true);
#endif

  // Read mesh
  ALE::Obj<Mesh> mesh;
  iohandler.read(&mesh);

  // Make sure meshIn matches data
  _checkVals(mesh, data);
} // _testRead

// ----------------------------------------------------------------------
// Test _orientCellsAscii with tet cells.
void
pylith::meshio::TestMeshIOLagrit::testOrientAsciiTet(void)
{ // testOrientAsciiTet
  // Expect vertices 1 and 2 to be swapped (not change to vertices 0 and 3)

  const int meshDim = 3;
  const int numCells = 2;
  const int numCorners = 4;
  const int cellsOrig[] = {
    0, 1, 2, 3,
    3, 4, 5, 6
  };
  const int cellsE[] = {
    0, 2, 1, 3,
    3, 5, 4, 6
  };
  
  int_array cells(cellsOrig, numCells*numCorners);
  MeshIOLagrit::_orientCellsAscii(&cells, numCells, numCorners, meshDim);

  const int size = numCells*numCorners;
  CPPUNIT_ASSERT_EQUAL(size, int(cells.size()));
  for (int i=0; i < size; ++i)
    CPPUNIT_ASSERT_EQUAL(cellsE[i], cells[i]);
} // testOrientAsciiTet

// ----------------------------------------------------------------------
// Test _orientCellsBinary with tet cells.
void
pylith::meshio::TestMeshIOLagrit::testOrientBinaryTet(void)
{ // testOrientBinaryTet
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
  MeshIOLagrit::_orientCellsBinary(&cells, numCells, numCorners, meshDim);

  const int size = numCells*numCorners;
  CPPUNIT_ASSERT_EQUAL(size, int(cells.size()));
  for (int i=0; i < size; ++i)
    CPPUNIT_ASSERT_EQUAL(cellsE[i], cells[i]);
} // testOrientBinaryTet


// End of file 
