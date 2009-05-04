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

#include "TestDirichletBCMultiTri3.hh" // Implementation of class methods

#include "data/DirichletDataMultiTri3.hh" // USES DirichletDataMultiTri3

// ----------------------------------------------------------------------
CPPUNIT_TEST_SUITE_REGISTRATION( pylith::bc::TestDirichletBCMultiTri3 );

// ----------------------------------------------------------------------
// Setup testing data.
void
pylith::bc::TestDirichletBCMultiTri3::setUp(void)
{ // setUp
  _data = new DirichletDataMultiTri3();
} // setUp


// End of file 
