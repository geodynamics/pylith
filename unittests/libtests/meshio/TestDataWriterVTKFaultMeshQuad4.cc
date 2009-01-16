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

#include "TestDataWriterVTKFaultMeshQuad4.hh" // Implementation of class methods

#include "data/DataWriterVTKDataFaultMeshQuad4.hh" // USES DataWriterVTKDataFaultMeshQuad4

// ----------------------------------------------------------------------
CPPUNIT_TEST_SUITE_REGISTRATION( pylith::meshio::TestDataWriterVTKFaultMeshQuad4 );

// ----------------------------------------------------------------------
// Setup testing data.
void
pylith::meshio::TestDataWriterVTKFaultMeshQuad4::setUp(void)
{ // setUp
  TestDataWriterVTKFaultMesh::setUp();
  _data = new DataWriterVTKDataFaultMeshQuad4;
  _flipFault = true;

  _initialize();
} // setUp


// End of file 
