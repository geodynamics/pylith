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

#include "TestDirichletBCTri3.hh" // Implementation of class methods

#include "data/DirichletDataTri3.hh" // USES DirichletDataTri3

// ----------------------------------------------------------------------
CPPUNIT_TEST_SUITE_REGISTRATION( pylith::bc::TestDirichletBCTri3 );

// ----------------------------------------------------------------------
// Setup testing data.
void
pylith::bc::TestDirichletBCTri3::setUp(void)
{ // setUp
  _data = new DirichletDataTri3();
} // setUp


// End of file 
