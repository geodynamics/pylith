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

#include "TestDataWriterVTKMeshTet4.hh" // Implementation of class methods

#include "data/DataWriterVTKDataMeshTet4.hh" // USES DataWriterVTKDataMeshTet4

// ----------------------------------------------------------------------
CPPUNIT_TEST_SUITE_REGISTRATION( pylith::meshio::TestDataWriterVTKMeshTet4 );

// ----------------------------------------------------------------------
// Setup testing data.
void
pylith::meshio::TestDataWriterVTKMeshTet4::setUp(void)
{ // setUp
  TestDataWriterVTKMesh::setUp();
  _data = new DataWriterVTKDataMeshTet4;
  _dataMesh = new DataWriterVTKDataMeshTet4;
  _initialize();
} // setUp


// End of file 
