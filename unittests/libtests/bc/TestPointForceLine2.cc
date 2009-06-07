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

#include "TestPointForceLine2.hh" // Implementation of class methods

#include "data/PointForceDataLine2.hh" // USES DirichletDataLine2

// ----------------------------------------------------------------------
CPPUNIT_TEST_SUITE_REGISTRATION( pylith::bc::TestPointForceLine2 );

// ----------------------------------------------------------------------
// Setup testing data.
void
pylith::bc::TestPointForceLine2::setUp(void)
{ // setUp
  _data = new PointForceDataLine2();
} // setUp


// End of file 
