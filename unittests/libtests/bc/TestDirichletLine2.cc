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

#include "TestDirichletLine2.hh" // Implementation of class methods

#include "data/DirichletDataLine2.hh" // USES DirichletDataLine2

// ----------------------------------------------------------------------
CPPUNIT_TEST_SUITE_REGISTRATION( pylith::bc::TestDirichletLine2 );

// ----------------------------------------------------------------------
// Setup testing data.
void
pylith::bc::TestDirichletLine2::setUp(void)
{ // setUp
  _data = new DirichletDataLine2();
} // setUp

// ----------------------------------------------------------------------
// Tear down testing data.
void
pylith::bc::TestDirichletLine2::tearDown(void)
{ // tearDown
  delete _data;
} // tearDown


// End of file 
