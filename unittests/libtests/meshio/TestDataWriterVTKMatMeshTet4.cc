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

#include "TestDataWriterVTKMatMeshTet4.hh" // Implementation of class methods

#include "data/DataWriterVTKDataMatMeshTet4.hh" // USES DataWriterVTKDataMeshTet4

// ----------------------------------------------------------------------
CPPUNIT_TEST_SUITE_REGISTRATION( pylith::meshio::TestDataWriterVTKMatMeshTet4 );

// ----------------------------------------------------------------------
// Setup testing data.
void
pylith::meshio::TestDataWriterVTKMatMeshTet4::setUp(void)
{ // setUp
  TestDataWriterVTKMesh::setUp();
  _data = new DataWriterVTKDataMatMeshTet4;
  _initialize();
} // setUp


// End of file 
