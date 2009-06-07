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

#include "TestPointForceTri3.hh" // Implementation of class methods

#include "data/PointForceDataTri3.hh" // USES DirichletDataTri3

// ----------------------------------------------------------------------
CPPUNIT_TEST_SUITE_REGISTRATION( pylith::bc::TestPointForceTri3 );

// ----------------------------------------------------------------------
// Setup testing data.
void
pylith::bc::TestPointForceTri3::setUp(void)
{ // setUp
  _data = new PointForceDataTri3();
} // setUp


// End of file 
