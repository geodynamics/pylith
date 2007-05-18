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

#include "TestDirichletTri3.hh" // Implementation of class methods

#include "data/DirichletDataTri3.hh" // USES DirichletDataTri3

// ----------------------------------------------------------------------
CPPUNIT_TEST_SUITE_REGISTRATION( pylith::bc::TestDirichletTri3 );

// ----------------------------------------------------------------------
// Setup testing data.
void
pylith::bc::TestDirichletTri3::setUp(void)
{ // setUp
  _data = new DirichletDataTri3();
} // setUp

// ----------------------------------------------------------------------
// Tear down testing data.
void
pylith::bc::TestDirichletTri3::tearDown(void)
{ // tearDown
  delete _data;
} // tearDown


// End of file 
