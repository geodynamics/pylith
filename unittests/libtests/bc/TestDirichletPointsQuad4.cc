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

#include "TestDirichletPointsQuad4.hh" // Implementation of class methods

#include "data/DirichletDataQuad4.hh" // USES DirichletDataQuad4

// ----------------------------------------------------------------------
CPPUNIT_TEST_SUITE_REGISTRATION( pylith::bc::TestDirichletPointsQuad4 );

// ----------------------------------------------------------------------
// Setup testing data.
void
pylith::bc::TestDirichletPointsQuad4::setUp(void)
{ // setUp
  _data = new DirichletDataQuad4();
} // setUp


// End of file 
