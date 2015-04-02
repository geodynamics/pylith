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
// Copyright (c) 2010-2015 University of California, Davis
//
// See COPYING for license information.
//
// ----------------------------------------------------------------------
//

#include <portinfo>

#include "TestPointForceCases.hh" // Implementation of class methods

// ----------------------------------------------------------------------
CPPUNIT_TEST_SUITE_REGISTRATION( pylith::bc::TestPointForceTri3 );

#include "data/PointForceDataTri3.hh" // USES DirichletDataTri3

// Setup testing data.
void
pylith::bc::TestPointForceTri3::setUp(void)
{ // setUp
  _data = new PointForceDataTri3();
} // setUp


// ----------------------------------------------------------------------
CPPUNIT_TEST_SUITE_REGISTRATION( pylith::bc::TestPointForceQuad4 );

#include "data/PointForceDataQuad4.hh" // USES DirichletDataQuad4

// Setup testing data.
void
pylith::bc::TestPointForceQuad4::setUp(void)
{ // setUp
  _data = new PointForceDataQuad4();
} // setUp


// ----------------------------------------------------------------------
CPPUNIT_TEST_SUITE_REGISTRATION( pylith::bc::TestPointForceTet4 );

#include "data/PointForceDataTet4.hh" // USES DirichletDataTet4

// Setup testing data.
void
pylith::bc::TestPointForceTet4::setUp(void)
{ // setUp
  _data = new PointForceDataTet4();
} // setUp


// ----------------------------------------------------------------------
CPPUNIT_TEST_SUITE_REGISTRATION( pylith::bc::TestPointForceHex8 );

#include "data/PointForceDataHex8.hh" // USES DirichletDataHex8

// Setup testing data.
void
pylith::bc::TestPointForceHex8::setUp(void)
{ // setUp
  _data = new PointForceDataHex8();
} // setUp


// End of file 
