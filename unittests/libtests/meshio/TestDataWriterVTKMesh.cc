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

#include "TestDataWriterVTKMesh.hh" // Implementation of class methods

#include "data/DataWriterVTKData.hh" // USES DataWriterVTKData
#include "data/DataWriterVTKDataMesh.hh" // USES DataWriterVTKDataMesh

#include "pylith/meshio/MeshIOAscii.hh" // USES MeshIOAscii

// ----------------------------------------------------------------------
// Setup testing data.
void
pylith::meshio::TestDataWriterVTKMesh::setUp(void)
{ // setUp
  TestDataWriterVTK::setUp();
  _dataMesh = 0;  
} // setUp

// ----------------------------------------------------------------------
// Tear down testing data.
void
pylith::meshio::TestDataWriterVTKMesh::tearDown(void)
{ // tearDown
  TestDataWriterVTK::tearDown();
  delete _dataMesh; _dataMesh = 0;
} // tearDown

// ----------------------------------------------------------------------
// Initialize mesh.
void
pylith::meshio::TestDataWriterVTKMesh::_initialize(void)
{ // _initialize
  CPPUNIT_ASSERT(0 != _dataMesh);

  MeshIOAscii iohandler;
  iohandler.filename(_dataMesh->meshFilename);
  iohandler.read(&_mesh);
  CPPUNIT_ASSERT(!_mesh.isNull());
  _mesh->getFactory()->clear();
} // _initialize


// End of file 
