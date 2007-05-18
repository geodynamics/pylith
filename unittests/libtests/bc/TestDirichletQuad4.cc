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

#include "TestDirichletQuad4.hh" // Implementation of class methods

#include "data/DirichletDataQuad4.hh" // USES DirichletDataQuad4

// ----------------------------------------------------------------------
CPPUNIT_TEST_SUITE_REGISTRATION( pylith::bc::TestDirichletQuad4 );

// ----------------------------------------------------------------------
// Setup testing data.
void
pylith::bc::TestDirichletQuad4::setUp(void)
{ // setUp
  _data = new DirichletDataQuad4();
} // setUp

// ----------------------------------------------------------------------
// Tear down testing data.
void
pylith::bc::TestDirichletQuad4::tearDown(void)
{ // tearDown
  delete _data;
} // tearDown


// End of file 
