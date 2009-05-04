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

#include "TestDirichletBCTet4.hh" // Implementation of class methods

#include "data/DirichletDataTet4.hh" // USES DirichletDataTet4

// ----------------------------------------------------------------------
CPPUNIT_TEST_SUITE_REGISTRATION( pylith::bc::TestDirichletBCTet4 );

// ----------------------------------------------------------------------
// Setup testing data.
void
pylith::bc::TestDirichletBCTet4::setUp(void)
{ // setUp
  _data = new DirichletDataTet4();
} // setUp


// End of file 
