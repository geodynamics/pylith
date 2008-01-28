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

#include "TestDirichletPointsTri3.hh" // Implementation of class methods

#include "data/DirichletPointsDataTri3.hh" // USES DirichletPointsDataTri3

// ----------------------------------------------------------------------
CPPUNIT_TEST_SUITE_REGISTRATION( pylith::bc::TestDirichletPointsTri3 );

// ----------------------------------------------------------------------
// Setup testing data.
void
pylith::bc::TestDirichletPointsTri3::setUp(void)
{ // setUp
  _data = new DirichletPointsDataTri3();
} // setUp


// End of file 
