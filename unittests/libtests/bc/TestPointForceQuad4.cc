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

#include "TestPointForceQuad4.hh" // Implementation of class methods

#include "data/PointForceDataQuad4.hh" // USES DirichletDataQuad4

// ----------------------------------------------------------------------
CPPUNIT_TEST_SUITE_REGISTRATION( pylith::bc::TestPointForceQuad4 );

// ----------------------------------------------------------------------
// Setup testing data.
void
pylith::bc::TestPointForceQuad4::setUp(void)
{ // setUp
  _data = new PointForceDataQuad4();
} // setUp


// End of file 
