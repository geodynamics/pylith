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

#include "TestDataWriterVTKSubMeshQuad4.hh" // Implementation of class methods

#include "data/DataWriterVTKDataSubMeshQuad4.hh" // USES DataWriterVTKDataMeshQuad4

// ----------------------------------------------------------------------
CPPUNIT_TEST_SUITE_REGISTRATION( pylith::meshio::TestDataWriterVTKSubMeshQuad4 );

// ----------------------------------------------------------------------
// Setup testing data.
void
pylith::meshio::TestDataWriterVTKSubMeshQuad4::setUp(void)
{ // setUp
  TestDataWriterVTKSubMesh::setUp();
  _data = new DataWriterVTKDataSubMeshQuad4;
  _initialize();
} // setUp


// End of file 
