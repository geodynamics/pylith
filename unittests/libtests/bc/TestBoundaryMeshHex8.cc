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
// Copyright (c) 2010-2011 University of California, Davis
//
// See COPYING for license information.
//
// ----------------------------------------------------------------------
//

#include <portinfo>

#include "TestBoundaryMeshHex8.hh" // Implementation of class methods

#include "data/BoundaryMeshDataHex8.hh" // USES BoundaryMeshDataHex8

// ----------------------------------------------------------------------
CPPUNIT_TEST_SUITE_REGISTRATION( pylith::bc::TestBoundaryMeshHex8 );

// ----------------------------------------------------------------------
// Setup testing data.
void
pylith::bc::TestBoundaryMeshHex8::setUp(void)
{ // setUp
  _data = new BoundaryMeshDataHex8();
  _flipFault = true;
} // setUp


// End of file 
