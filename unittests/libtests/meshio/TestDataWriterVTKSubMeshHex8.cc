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

#include "TestDataWriterVTKSubMeshHex8.hh" // Implementation of class methods

#include "data/DataWriterVTKDataSubMeshHex8.hh" // USES DataWriterVTKDataMeshHex8

// ----------------------------------------------------------------------
CPPUNIT_TEST_SUITE_REGISTRATION( pylith::meshio::TestDataWriterVTKSubMeshHex8 );

// ----------------------------------------------------------------------
// Setup testing data.
void
pylith::meshio::TestDataWriterVTKSubMeshHex8::setUp(void)
{ // setUp
  TestDataWriterVTKSubMesh::setUp();
  _data = new DataWriterVTKDataSubMeshHex8;
  _flipFault = true;

  _initialize();
} // setUp


// End of file 
