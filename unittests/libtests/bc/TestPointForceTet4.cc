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

#include "TestPointForceTet4.hh" // Implementation of class methods

#include "data/PointForceDataTet4.hh" // USES DirichletDataTet4

// ----------------------------------------------------------------------
CPPUNIT_TEST_SUITE_REGISTRATION( pylith::bc::TestPointForceTet4 );

// ----------------------------------------------------------------------
// Setup testing data.
void
pylith::bc::TestPointForceTet4::setUp(void)
{ // setUp
  _data = new PointForceDataTet4();
} // setUp


// End of file 
