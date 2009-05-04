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

#include "TestDataWriterVTKSubMeshLine2.hh" // Implementation of class methods

#include "data/DataWriterVTKDataSubMeshLine2.hh" // USES DataWriterVTKDataMeshLine2

// ----------------------------------------------------------------------
CPPUNIT_TEST_SUITE_REGISTRATION( pylith::meshio::TestDataWriterVTKSubMeshLine2 );

// ----------------------------------------------------------------------
// Setup testing data.
void
pylith::meshio::TestDataWriterVTKSubMeshLine2::setUp(void)
{ // setUp
  TestDataWriterVTKSubMesh::setUp();
  _data = new DataWriterVTKDataSubMeshLine2;
  _initialize();
} // setUp


// End of file 
