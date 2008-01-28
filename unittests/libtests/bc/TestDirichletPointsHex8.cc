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

#include "TestDirichletPointsHex8.hh" // Implementation of class methods

#include "data/DirichletPointsDataHex8.hh" // USES DirichletPointsDataHex8

// ----------------------------------------------------------------------
CPPUNIT_TEST_SUITE_REGISTRATION( pylith::bc::TestDirichletPointsHex8 );

// ----------------------------------------------------------------------
// Setup testing data.
void
pylith::bc::TestDirichletPointsHex8::setUp(void)
{ // setUp
  _data = new DirichletPointsDataHex8();
} // setUp


// End of file 
