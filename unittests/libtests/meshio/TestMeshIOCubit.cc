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


// End of file 
