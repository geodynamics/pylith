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

#include "TestDirichletPointsLine2.hh" // Implementation of class methods

#include "data/DirichletDataLine2.hh" // USES DirichletDataLine2

// ----------------------------------------------------------------------
CPPUNIT_TEST_SUITE_REGISTRATION( pylith::bc::TestDirichletPointsLine2 );

// ----------------------------------------------------------------------
// Setup testing data.
void
pylith::bc::TestDirichletPointsLine2::setUp(void)
{ // setUp
  _data = new DirichletDataLine2();
} // setUp


// End of file 
