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

#include "TestGeometryLine1D.hh" // Implementation of class methods

#include "pylith/feassemble/GeometryLine1D.hh"
#include "pylith/feassemble/GeometryPoint1D.hh"
#include "pylith/feassemble/GeometryPoint2D.hh"

#include "data/GeomDataLine1D.hh"

// ----------------------------------------------------------------------
CPPUNIT_TEST_SUITE_REGISTRATION( pylith::feassemble::TestGeometryLine1D );

// ----------------------------------------------------------------------
// Setup test data.
void
pylith::feassemble::TestGeometryLine1D::setUp(void)
{ // setUp
  _object = new GeometryLine1D();
  _data = new GeomDataLine1D();
} // setUp

// ----------------------------------------------------------------------
// Test constructor
void
pylith::feassemble::TestGeometryLine1D::testConstructor(void)
{ // testConstructor
  GeometryLine1D geometry;
} // testConstructor

// ----------------------------------------------------------------------
// Test geometryLowerDim().
void
pylith::feassemble::TestGeometryLine1D::testGeomLowerDim(void)
{ // testGeomLowerDim
  GeometryLine1D geometry;
  CellGeometry* geometryLD = geometry.geometryLowerDim();
  GeometryPoint1D* geometryPt = dynamic_cast<GeometryPoint1D*>(geometryLD);
  CPPUNIT_ASSERT(0 != geometryPt);
  GeometryPoint2D* geometryPt2 = dynamic_cast<GeometryPoint2D*>(geometryLD);
  CPPUNIT_ASSERT(0 == geometryPt2);
  delete geometryLD; geometryLD = 0;
} // testGeomLowerDim


// End of file 
