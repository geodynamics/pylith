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

#include "TestBoundaryMeshTri3.hh" // Implementation of class methods

#include "data/BoundaryMeshDataTri3.hh" // USES BoundaryMeshDataTri3

#include "pylith/faults/CohesiveTopology.hh"

// ----------------------------------------------------------------------
CPPUNIT_TEST_SUITE_REGISTRATION( pylith::bc::TestBoundaryMeshTri3 );

// ----------------------------------------------------------------------
// Setup testing data.
void
pylith::bc::TestBoundaryMeshTri3::setUp(void)
{ // setUp
  _data = new BoundaryMeshDataTri3();
  _flipFault = true;
} // setUp


// End of file 
