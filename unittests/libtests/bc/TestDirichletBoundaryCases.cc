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

#include "TestDirichletBoundaryCases.hh" // Implementation of cases

// ----------------------------------------------------------------------
CPPUNIT_TEST_SUITE_REGISTRATION( pylith::bc::TestDirichletBoundaryTri3 );

#include "data/DirichletDataTri3.hh" // USES DirichletDataTri3

// Setup testing data.
void
pylith::bc::TestDirichletBoundaryTri3::setUp(void)
{ // setUp
  _data = new DirichletDataTri3();
} // setUp


// ----------------------------------------------------------------------
CPPUNIT_TEST_SUITE_REGISTRATION( pylith::bc::TestDirichletBoundaryQuad4 );

#include "data/DirichletDataQuad4.hh" // USES DirichletDataQuad4

// Setup testing data.
void
pylith::bc::TestDirichletBoundaryQuad4::setUp(void)
{ // setUp
  _data = new DirichletDataQuad4();
} // setUp


// ----------------------------------------------------------------------
CPPUNIT_TEST_SUITE_REGISTRATION( pylith::bc::TestDirichletBoundaryTet4 );

#include "data/DirichletDataTet4.hh" // USES DirichletDataTet4

// Setup testing data.
void
pylith::bc::TestDirichletBoundaryTet4::setUp(void)
{ // setUp
  _data = new DirichletDataTet4();
} // setUp


// ----------------------------------------------------------------------
CPPUNIT_TEST_SUITE_REGISTRATION( pylith::bc::TestDirichletBoundaryHex8 );

#include "data/DirichletDataHex8.hh" // USES DirichletDataHex8

// Setup testing data.
void
pylith::bc::TestDirichletBoundaryHex8::setUp(void)
{ // setUp
  _data = new DirichletDataHex8();
} // setUp


// End of file 
