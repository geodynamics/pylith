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

#include "TestBoundaryMeshCases.hh" // Implementation of cases

#include "pylith/faults/CohesiveTopology.hh"

// ----------------------------------------------------------------------
CPPUNIT_TEST_SUITE_REGISTRATION( pylith::bc::TestBoundaryMeshTri3 );

#include "data/BoundaryMeshDataTri3.hh" // USES BoundaryMeshDataTri3

// Setup testing data.
void
pylith::bc::TestBoundaryMeshTri3::setUp(void)
{ // setUp
  _data = new BoundaryMeshDataTri3();
} // setUp


// ----------------------------------------------------------------------
CPPUNIT_TEST_SUITE_REGISTRATION( pylith::bc::TestBoundaryMeshQuad4 );

#include "data/BoundaryMeshDataQuad4.hh" // USES BoundaryMeshDataQuad4

// Setup testing data.
void
pylith::bc::TestBoundaryMeshQuad4::setUp(void)
{ // setUp
  _data = new BoundaryMeshDataQuad4();
} // setUp


// ----------------------------------------------------------------------
CPPUNIT_TEST_SUITE_REGISTRATION( pylith::bc::TestBoundaryMeshTet4 );

#include "data/BoundaryMeshDataTet4.hh" // USES BoundaryMeshDataTet4

// Setup testing data.
void
pylith::bc::TestBoundaryMeshTet4::setUp(void)
{ // setUp
  _data = new BoundaryMeshDataTet4();
} // setUp


// ----------------------------------------------------------------------
CPPUNIT_TEST_SUITE_REGISTRATION( pylith::bc::TestBoundaryMeshHex8 );

#include "data/BoundaryMeshDataHex8.hh" // USES BoundaryMeshDataHex8

// Setup testing data.
void
pylith::bc::TestBoundaryMeshHex8::setUp(void)
{ // setUp
  _data = new BoundaryMeshDataHex8();
} // setUp


// End of file 
