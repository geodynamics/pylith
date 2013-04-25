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
// Copyright (c) 2010-2013 University of California, Davis
//
// See COPYING for license information.
//
// ----------------------------------------------------------------------
//

#include <portinfo>

#include "TestDirichletBCCases.hh" // Implementation of cases

// ----------------------------------------------------------------------
CPPUNIT_TEST_SUITE_REGISTRATION( pylith::bc::TestDirichletBCLine2 );

#include "data/DirichletDataLine2.hh" // USES DirichletDataLine2

// Setup testing data.
void
pylith::bc::TestDirichletBCLine2::setUp(void)
{ // setUp
  _data = new DirichletDataLine2();
} // setUp


// ----------------------------------------------------------------------
CPPUNIT_TEST_SUITE_REGISTRATION( pylith::bc::TestDirichletBCLine2b );

#include "data/DirichletDataLine2b.hh" // USES DirichletDataLine2b

// Setup testing data.
void
pylith::bc::TestDirichletBCLine2b::setUp(void)
{ // setUp
  _data = new DirichletDataLine2b();
} // setUp


// ----------------------------------------------------------------------
CPPUNIT_TEST_SUITE_REGISTRATION( pylith::bc::TestDirichletBCTri3 );

#include "data/DirichletDataTri3.hh" // USES DirichletDataTri3

// Setup testing data.
void
pylith::bc::TestDirichletBCTri3::setUp(void)
{ // setUp
  _data = new DirichletDataTri3();
} // setUp


// ----------------------------------------------------------------------
CPPUNIT_TEST_SUITE_REGISTRATION( pylith::bc::TestDirichletBCQuad4 );

#include "data/DirichletDataQuad4.hh" // USES DirichletDataQuad4

// Setup testing data.
void
pylith::bc::TestDirichletBCQuad4::setUp(void)
{ // setUp
  _data = new DirichletDataQuad4();
} // setUp


// ----------------------------------------------------------------------
CPPUNIT_TEST_SUITE_REGISTRATION( pylith::bc::TestDirichletBCTet4 );

#include "data/DirichletDataTet4.hh" // USES DirichletDataTet4

// Setup testing data.
void
pylith::bc::TestDirichletBCTet4::setUp(void)
{ // setUp
  _data = new DirichletDataTet4();
} // setUp


// ----------------------------------------------------------------------
CPPUNIT_TEST_SUITE_REGISTRATION( pylith::bc::TestDirichletBCHex8 );

#include "data/DirichletDataHex8.hh" // USES DirichletDataHex8

// Setup testing data.
void
pylith::bc::TestDirichletBCHex8::setUp(void)
{ // setUp
  _data = new DirichletDataHex8();
} // setUp


// End of file 
