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

#include "TestDirichletBCMultiTet4.hh" // Implementation of class methods

#include "data/DirichletDataMultiTet4.hh" // USES DirichletDataMultiTet4

// ----------------------------------------------------------------------
CPPUNIT_TEST_SUITE_REGISTRATION( pylith::bc::TestDirichletBCMultiTet4 );

// ----------------------------------------------------------------------
// Setup testing data.
void
pylith::bc::TestDirichletBCMultiTet4::setUp(void)
{ // setUp
  _data = new DirichletDataMultiTet4();
} // setUp


// End of file 
