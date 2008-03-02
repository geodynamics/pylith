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

#include "TestBoundaryMeshHex8.hh" // Implementation of class methods

#include "data/BoundaryMeshDataHex8.hh" // USES BoundaryMeshDataHex8

// ----------------------------------------------------------------------
CPPUNIT_TEST_SUITE_REGISTRATION( pylith::bc::TestBoundaryMeshHex8 );

// ----------------------------------------------------------------------
// Setup testing data.
void
pylith::bc::TestBoundaryMeshHex8::setUp(void)
{ // setUp
  _data = new BoundaryMeshDataHex8();
} // setUp


// End of file 
