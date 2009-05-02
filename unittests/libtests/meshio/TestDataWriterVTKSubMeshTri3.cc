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

#include "TestDataWriterVTKSubMeshTri3.hh" // Implementation of class methods

#include "data/DataWriterVTKDataSubMeshTri3.hh" // USES DataWriterVTKDataMeshTri3

// ----------------------------------------------------------------------
CPPUNIT_TEST_SUITE_REGISTRATION( pylith::meshio::TestDataWriterVTKSubMeshTri3 );

// ----------------------------------------------------------------------
// Setup testing data.
void
pylith::meshio::TestDataWriterVTKSubMeshTri3::setUp(void)
{ // setUp
  TestDataWriterVTKSubMesh::setUp();
  _data = new DataWriterVTKDataSubMeshTri3;
  _flipFault = true;

  _initialize();
} // setUp


// End of file 
