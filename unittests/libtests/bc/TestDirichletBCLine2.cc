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

#include "TestDirichletBCLine2.hh" // Implementation of class methods

#include "data/DirichletDataLine2.hh" // USES DirichletDataLine2

// ----------------------------------------------------------------------
CPPUNIT_TEST_SUITE_REGISTRATION( pylith::bc::TestDirichletBCLine2 );

// ----------------------------------------------------------------------
// Setup testing data.
void
pylith::bc::TestDirichletBCLine2::setUp(void)
{ // setUp
  _data = new DirichletDataLine2();
} // setUp


// End of file 
