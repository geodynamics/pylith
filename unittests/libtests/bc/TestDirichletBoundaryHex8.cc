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

#include "TestDirichletBoundaryHex8.hh" // Implementation of class methods

#include "data/DirichletDataHex8.hh" // USES DirichletDataHex8

// ----------------------------------------------------------------------
CPPUNIT_TEST_SUITE_REGISTRATION( pylith::bc::TestDirichletBoundaryHex8 );

// ----------------------------------------------------------------------
// Setup testing data.
void
pylith::bc::TestDirichletBoundaryHex8::setUp(void)
{ // setUp
  _data = new DirichletDataHex8();
} // setUp


// End of file 
