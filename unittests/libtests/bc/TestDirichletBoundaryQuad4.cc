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

#include "TestDirichletBoundaryQuad4.hh" // Implementation of class methods

#include "data/DirichletDataQuad4.hh" // USES DirichletDataQuad4

// ----------------------------------------------------------------------
CPPUNIT_TEST_SUITE_REGISTRATION( pylith::bc::TestDirichletBoundaryQuad4 );

// ----------------------------------------------------------------------
// Setup testing data.
void
pylith::bc::TestDirichletBoundaryQuad4::setUp(void)
{ // setUp
  _data = new DirichletDataQuad4();
} // setUp


// End of file 
