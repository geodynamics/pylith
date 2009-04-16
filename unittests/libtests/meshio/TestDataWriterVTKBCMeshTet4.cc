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

#include "TestDataWriterVTKBCMeshTet4.hh" // Implementation of class methods

#include "data/DataWriterVTKDataBCMeshTet4.hh" // USES DataWriterVTKDataBCMeshTet4

// ----------------------------------------------------------------------
CPPUNIT_TEST_SUITE_REGISTRATION( pylith::meshio::TestDataWriterVTKBCMeshTet4 );

// ----------------------------------------------------------------------
// Setup testing data.
void
pylith::meshio::TestDataWriterVTKBCMeshTet4::setUp(void)
{ // setUp
  TestDataWriterVTKBCMesh::setUp();
  _data = new DataWriterVTKDataBCMeshTet4;
  _initialize();
} // setUp


// End of file 
