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

#include "TestBoundaryMeshQuad4.hh" // Implementation of class methods

#include "data/BoundaryMeshDataQuad4.hh" // USES BoundaryMeshDataQuad4

// ----------------------------------------------------------------------
CPPUNIT_TEST_SUITE_REGISTRATION( pylith::bc::TestBoundaryMeshQuad4 );

// ----------------------------------------------------------------------
// Setup testing data.
void
pylith::bc::TestBoundaryMeshQuad4::setUp(void)
{ // setUp
  _data = new BoundaryMeshDataQuad4();
  _flipFault = true;
} // setUp


// End of file 
