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

#include "TestPointForceHex8.hh" // Implementation of class methods

#include "data/PointForceDataHex8.hh" // USES DirichletDataHex8

// ----------------------------------------------------------------------
CPPUNIT_TEST_SUITE_REGISTRATION( pylith::bc::TestPointForceHex8 );

// ----------------------------------------------------------------------
// Setup testing data.
void
pylith::bc::TestPointForceHex8::setUp(void)
{ // setUp
  _data = new PointForceDataHex8();
} // setUp


// End of file 
