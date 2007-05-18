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

#include "TestDirichletHex8.hh" // Implementation of class methods

#include "data/DirichletDataHex8.hh" // USES DirichletDataHex8

// ----------------------------------------------------------------------
CPPUNIT_TEST_SUITE_REGISTRATION( pylith::bc::TestDirichletHex8 );

// ----------------------------------------------------------------------
// Setup testing data.
void
pylith::bc::TestDirichletHex8::setUp(void)
{ // setUp
  _data = new DirichletDataHex8();
} // setUp

// ----------------------------------------------------------------------
// Tear down testing data.
void
pylith::bc::TestDirichletHex8::tearDown(void)
{ // tearDown
  delete _data;
} // tearDown


// End of file 
