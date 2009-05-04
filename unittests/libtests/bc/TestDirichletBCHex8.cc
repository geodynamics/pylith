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

#include "TestDirichletBCHex8.hh" // Implementation of class methods

#include "data/DirichletDataHex8.hh" // USES DirichletDataHex8

// ----------------------------------------------------------------------
CPPUNIT_TEST_SUITE_REGISTRATION( pylith::bc::TestDirichletBCHex8 );

// ----------------------------------------------------------------------
// Setup testing data.
void
pylith::bc::TestDirichletBCHex8::setUp(void)
{ // setUp
  _data = new DirichletDataHex8();
} // setUp


// End of file 
