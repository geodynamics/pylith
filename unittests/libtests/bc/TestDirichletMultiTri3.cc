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

#include "TestDirichletMultiTri3.hh" // Implementation of class methods

#include "data/DirichletDataMultiTri3.hh" // USES DirichletDataMultiTri3

// ----------------------------------------------------------------------
CPPUNIT_TEST_SUITE_REGISTRATION( pylith::bc::TestDirichletMultiTri3 );

// ----------------------------------------------------------------------
// Setup testing data.
void
pylith::bc::TestDirichletMultiTri3::setUp(void)
{ // setUp
  _data = new DirichletDataMultiTri3();
} // setUp


// End of file 
