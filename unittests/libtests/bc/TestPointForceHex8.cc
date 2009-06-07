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

#include "TestPointForceHex8.hh" // Implementation of class methods

#include "data/PointForceDataHex8.hh" // USES DirichletDataHex8

// ----------------------------------------------------------------------
CPPUNIT_TEST_SUITE_REGISTRATION( pylith::bc::TestPointForceHex8 );

// ----------------------------------------------------------------------
// Setup testing data.
void
pylith::bc::TestPointForceHex8::setUp(void)
{ // setUp
  _data = new PointForceDataHex8();
} // setUp


// End of file 
