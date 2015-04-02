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

#include "TestDirichletBCMultiCases.hh" // Implementation of cases

// ----------------------------------------------------------------------
CPPUNIT_TEST_SUITE_REGISTRATION( pylith::bc::TestDirichletBCMultiTri3 );

#include "data/DirichletDataMultiTri3.hh" // USES DirichletDataMultiTri3

// Setup testing data.
void
pylith::bc::TestDirichletBCMultiTri3::setUp(void)
{ // setUp
  _data = new DirichletDataMultiTri3();
} // setUp


// ----------------------------------------------------------------------
CPPUNIT_TEST_SUITE_REGISTRATION( pylith::bc::TestDirichletBCMultiTet4 );

#include "data/DirichletDataMultiTet4.hh" // USES DirichletDataMultiTet4

// Setup testing data.
void
pylith::bc::TestDirichletBCMultiTet4::setUp(void)
{ // setUp
  _data = new DirichletDataMultiTet4();
} // setUp


// End of file 
