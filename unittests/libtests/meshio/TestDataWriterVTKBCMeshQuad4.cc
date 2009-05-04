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

#include "TestDataWriterVTKBCMeshQuad4.hh" // Implementation of class methods

#include "data/DataWriterVTKDataBCMeshQuad4.hh" // USES DataWriterVTKDataMeshQuad4

// ----------------------------------------------------------------------
CPPUNIT_TEST_SUITE_REGISTRATION( pylith::meshio::TestDataWriterVTKBCMeshQuad4 );

// ----------------------------------------------------------------------
// Setup testing data.
void
pylith::meshio::TestDataWriterVTKBCMeshQuad4::setUp(void)
{ // setUp
  TestDataWriterVTKBCMesh::setUp();
  _data = new DataWriterVTKDataBCMeshQuad4;
  _initialize();
} // setUp


// End of file 
