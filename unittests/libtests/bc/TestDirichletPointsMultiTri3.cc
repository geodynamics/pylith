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

#include "TestDirichletPointsMultiTri3.hh" // Implementation of class methods

#include "data/DirichletPointsDataMultiTri3.hh" // USES DirichletPointsDataMultiTri3

// ----------------------------------------------------------------------
CPPUNIT_TEST_SUITE_REGISTRATION( pylith::bc::TestDirichletPointsMultiTri3 );

// ----------------------------------------------------------------------
// Setup testing data.
void
pylith::bc::TestDirichletPointsMultiTri3::setUp(void)
{ // setUp
  _data = new DirichletPointsDataMultiTri3();
} // setUp


// End of file 
