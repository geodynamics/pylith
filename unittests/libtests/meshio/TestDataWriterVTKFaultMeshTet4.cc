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

#include "TestDataWriterVTKFaultMeshTet4.hh" // Implementation of class methods

#include "data/DataWriterVTKDataFaultMeshTet4.hh" // USES DataWriterVTKDataFaultMeshTet4

// ----------------------------------------------------------------------
CPPUNIT_TEST_SUITE_REGISTRATION( pylith::meshio::TestDataWriterVTKFaultMeshTet4 );

// ----------------------------------------------------------------------
// Setup testing data.
void
pylith::meshio::TestDataWriterVTKFaultMeshTet4::setUp(void)
{ // setUp
  TestDataWriterVTKFaultMesh::setUp();
  _data = new DataWriterVTKDataFaultMeshTet4;
  _dataMesh = new DataWriterVTKDataFaultMeshTet4;
  _initialize();
} // setUp


// End of file 
