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

#include "TestDataWriterVTKBCMeshHex8.hh" // Implementation of class methods

#include "data/DataWriterVTKDataBCMeshHex8.hh" // USES DataWriterVTKDataBCMeshHex8

// ----------------------------------------------------------------------
CPPUNIT_TEST_SUITE_REGISTRATION( pylith::meshio::TestDataWriterVTKBCMeshHex8 );

// ----------------------------------------------------------------------
// Setup testing data.
void
pylith::meshio::TestDataWriterVTKBCMeshHex8::setUp(void)
{ // setUp
  TestDataWriterVTKBCMesh::setUp();
  _data = new DataWriterVTKDataBCMeshHex8;
  _flipFault = true;
  _initialize();
} // setUp


// End of file 
