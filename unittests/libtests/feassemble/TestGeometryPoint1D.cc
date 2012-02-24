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

#include "TestGeometryPoint1D.hh" // Implementation of class methods

#include "pylith/feassemble/GeometryPoint1D.hh"

#include "pylith/utils/array.hh" // USES scalar_array

#include "data/GeomDataPoint1D.hh"

// ----------------------------------------------------------------------
CPPUNIT_TEST_SUITE_REGISTRATION( pylith::feassemble::TestGeometryPoint1D );

// ----------------------------------------------------------------------
// Setup test data.
void
pylith::feassemble::TestGeometryPoint1D::setUp(void)
{ // setUp
  _object = new GeometryPoint1D();
  _data = new GeomDataPoint1D();
} // setUp

// ----------------------------------------------------------------------
// Test constructor
void
pylith::feassemble::TestGeometryPoint1D::testConstructor(void)
{ // testConstructor
  GeometryPoint1D geometry;
} // testConstructor

// ----------------------------------------------------------------------
// Test geometryLowerDim().
void
pylith::feassemble::TestGeometryPoint1D::testGeomLowerDim(void)
{ // testGeomLowerDim
  GeometryPoint1D geometry;
  CellGeometry* geometryLD = geometry.geometryLowerDim();
  CPPUNIT_ASSERT(0 == geometryLD);
  delete geometryLD; geometryLD = 0;
} // testGeomLowerDim

// ----------------------------------------------------------------------
// Test jacobian().
void
pylith::feassemble::TestGeometryPoint1D::testJacobian(void)
{ // testJacobian
  GeometryPoint1D geometry;
  GeomDataPoint1D data;

  const int cellDim = data.cellDim;
  const int spaceDim = data.spaceDim;
  const int numCorners = data.numCorners;

  const int numLocs = data.numLocs;

  CPPUNIT_ASSERT_EQUAL(cellDim, geometry.cellDim());
  CPPUNIT_ASSERT_EQUAL(spaceDim, geometry.spaceDim());
  CPPUNIT_ASSERT_EQUAL(numCorners, geometry.numCorners());

  scalar_array jacobian(1);
  PylithScalar det = 0.0;
  for (int iLoc=0; iLoc < numLocs; ++iLoc) {
    scalar_array vertices(data.vertices, numCorners*spaceDim);
    scalar_array location(&data.locations[iLoc], 1);

    geometry.jacobian(&jacobian, &det, vertices, location);
    CPPUNIT_ASSERT_EQUAL(PylithScalar(1.0), jacobian[0]);
    CPPUNIT_ASSERT_EQUAL(PylithScalar(1.0), det);
  } //for
} // testJacobian


// End of file 
