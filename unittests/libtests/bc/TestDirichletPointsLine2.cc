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

#include "data/DirichletPointsDataLine2.hh" // USES DirichletPointsDataLine2

// ----------------------------------------------------------------------
CPPUNIT_TEST_SUITE_REGISTRATION( pylith::bc::TestDirichletPointsLine2 );

// ----------------------------------------------------------------------
// Setup testing data.
void
pylith::bc::TestDirichletPointsLine2::setUp(void)
{ // setUp
  _data = new DirichletPointsDataLine2();
} // setUp


// End of file 
