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

#include "TestDirichletPointsMultiTet4.hh" // Implementation of class methods

#include "data/DirichletDataMultiTet4.hh" // USES DirichletDataMultiTet4

// ----------------------------------------------------------------------
CPPUNIT_TEST_SUITE_REGISTRATION( pylith::bc::TestDirichletPointsMultiTet4 );

// ----------------------------------------------------------------------
// Setup testing data.
void
pylith::bc::TestDirichletPointsMultiTet4::setUp(void)
{ // setUp
  _data = new DirichletDataMultiTet4();
} // setUp


// End of file 
