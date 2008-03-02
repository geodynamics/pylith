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

#include "TestBoundaryMeshTet4.hh" // Implementation of class methods

#include "data/BoundaryMeshDataTet4.hh" // USES BoundaryMeshDataTet4

// ----------------------------------------------------------------------
CPPUNIT_TEST_SUITE_REGISTRATION( pylith::bc::TestBoundaryMeshTet4 );

// ----------------------------------------------------------------------
// Setup testing data.
void
pylith::bc::TestBoundaryMeshTet4::setUp(void)
{ // setUp
  _data = new BoundaryMeshDataTet4();
} // setUp


// End of file 
