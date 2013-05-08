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

#include "TestGeometryPoint2D.hh" // Implementation of class methods

#include "pylith/feassemble/GeometryPoint2D.hh"

#include "pylith/utils/array.hh" // USES scalar_array

#include "data/GeomDataPoint2D.hh"

#include "pylith/utils/error.h" // USES PYLITH_METHOD_BEGIN/END

// ----------------------------------------------------------------------
CPPUNIT_TEST_SUITE_REGISTRATION( pylith::feassemble::TestGeometryPoint2D );

// ----------------------------------------------------------------------
// Setup test data.
void
pylith::feassemble::TestGeometryPoint2D::setUp(void)
{ // setUp
  PYLITH_METHOD_BEGIN;
  
  _object = new GeometryPoint2D();
  _data = new GeomDataPoint2D();

  PYLITH_METHOD_END;
} // setUp

// ----------------------------------------------------------------------
// Test constructor
void
pylith::feassemble::TestGeometryPoint2D::testConstructor(void)
{ // testConstructor
  PYLITH_METHOD_BEGIN;
  
  GeometryPoint2D geometry;

  PYLITH_METHOD_END;
} // testConstructor

// ----------------------------------------------------------------------
// Test geometryLowerDim().
void
pylith::feassemble::TestGeometryPoint2D::testGeomLowerDim(void)
{ // testGeomLowerDim
  PYLITH_METHOD_BEGIN;
  
  GeometryPoint2D geometry;
  CellGeometry* geometryLD = geometry.geometryLowerDim();
  CPPUNIT_ASSERT(0 == geometryLD);
  delete geometryLD; geometryLD = 0;

  PYLITH_METHOD_END;
} // testGeomLowerDim

// ----------------------------------------------------------------------
// Test jacobian().
void
pylith::feassemble::TestGeometryPoint2D::testJacobian(void)
{ // testJacobian
  PYLITH_METHOD_BEGIN;
  
  GeometryPoint2D geometry;
  GeomDataPoint2D data;

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
    geometry.jacobian(&jacobian, &det, data.vertices, numCorners, spaceDim, &data.locations[iLoc], 1);
    CPPUNIT_ASSERT_EQUAL(PylithScalar(1.0), jacobian[0]);
    CPPUNIT_ASSERT_EQUAL(PylithScalar(1.0), det);
  } //for

  PYLITH_METHOD_END;
} // testJacobian


// End of file 
