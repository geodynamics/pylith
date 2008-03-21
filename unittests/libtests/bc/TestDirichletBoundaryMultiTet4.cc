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

#include "TestDirichletBoundaryMultiTet4.hh" // Implementation of class methods

#include "data/DirichletDataMultiTet4.hh" // USES DirichletDataMultiTet4

// ----------------------------------------------------------------------
CPPUNIT_TEST_SUITE_REGISTRATION( pylith::bc::TestDirichletBoundaryMultiTet4 );

// ----------------------------------------------------------------------
// Setup testing data.
void
pylith::bc::TestDirichletBoundaryMultiTet4::setUp(void)
{ // setUp
  _data = new DirichletDataMultiTet4();
} // setUp


// End of file 
