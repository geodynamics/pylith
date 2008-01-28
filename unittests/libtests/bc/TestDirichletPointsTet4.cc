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

#include "TestDirichletPointsTet4.hh" // Implementation of class methods

#include "data/DirichletPointsDataTet4.hh" // USES DirichletPointsDataTet4

// ----------------------------------------------------------------------
CPPUNIT_TEST_SUITE_REGISTRATION( pylith::bc::TestDirichletPointsTet4 );

// ----------------------------------------------------------------------
// Setup testing data.
void
pylith::bc::TestDirichletPointsTet4::setUp(void)
{ // setUp
  _data = new DirichletPointsDataTet4();
} // setUp


// End of file 
