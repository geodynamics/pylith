// -*- C++ -*-
//
// ----------------------------------------------------------------------
//
// Brad T. Aagaard, U.S. Geological Survey
// Charles A. Williams, GNS Science
// Matthew G. Knepley, University of Chicago
//
// This code was developed as part of the Computational Infrastructure
// for Geodynamics (http://geodynamics.org).
//
// Copyright (c) 2010-2012 University of California, Davis
//
// See COPYING for license information.
//
// ----------------------------------------------------------------------
//

#include <portinfo>

#include "TestPointForceQuad4.hh" // Implementation of class methods

#include "data/PointForceDataQuad4.hh" // USES DirichletDataQuad4

// ----------------------------------------------------------------------
CPPUNIT_TEST_SUITE_REGISTRATION( pylith::bc::TestPointForceQuad4 );

// ----------------------------------------------------------------------
// Setup testing data.
void
pylith::bc::TestPointForceQuad4::setUp(void)
{ // setUp
  _data = new PointForceDataQuad4();
} // setUp


// End of file 
