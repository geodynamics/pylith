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

#include "TestBoundaryMeshTri3.hh" // Implementation of class methods

#include "data/BoundaryMeshDataTri3.hh" // USES BoundaryMeshDataTri3

#include "pylith/faults/CohesiveTopology.hh"

// ----------------------------------------------------------------------
CPPUNIT_TEST_SUITE_REGISTRATION( pylith::bc::TestBoundaryMeshTri3 );

// ----------------------------------------------------------------------
// Setup testing data.
void
pylith::bc::TestBoundaryMeshTri3::setUp(void)
{ // setUp
  _data = new BoundaryMeshDataTri3();
  _flipFault = true;
} // setUp


// End of file 
