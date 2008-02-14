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

#include "TestDirichletPointsLine2b.hh" // Implementation of class methods

#include "data/DirichletDataLine2b.hh" // USES DirichletDataLine2b

// ----------------------------------------------------------------------
CPPUNIT_TEST_SUITE_REGISTRATION( pylith::bc::TestDirichletPointsLine2b );

// ----------------------------------------------------------------------
// Setup testing data.
void
pylith::bc::TestDirichletPointsLine2b::setUp(void)
{ // setUp
  _data = new DirichletDataLine2b();
} // setUp


// End of file 
