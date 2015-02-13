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

#include "TestGeometryLine2D.hh" // Implementation of class methods

#include "pylith/feassemble/GeometryLine2D.hh"

#include "data/GeomDataLine2D.hh"

#include "pylith/utils/error.h" // USES PYLITH_METHOD_BEGIN/END

// ----------------------------------------------------------------------
CPPUNIT_TEST_SUITE_REGISTRATION( pylith::feassemble::TestGeometryLine2D );

// ----------------------------------------------------------------------
// Setup test data.
void
pylith::feassemble::TestGeometryLine2D::setUp(void)
{ // setUp
  PYLITH_METHOD_BEGIN;

  _object = new GeometryLine2D();
  _data = new GeomDataLine2D();

  PYLITH_METHOD_END;
} // setUp

// ----------------------------------------------------------------------
// Test constructor
void
pylith::feassemble::TestGeometryLine2D::testConstructor(void)
{ // testConstructor
  PYLITH_METHOD_BEGIN;

  GeometryLine2D geometry;

  PYLITH_METHOD_END;
} // testConstructor

// ----------------------------------------------------------------------
// Test geometryLowerDim().
void
pylith::feassemble::TestGeometryLine2D::testGeomLowerDim(void)
{ // testGeomLowerDim
  PYLITH_METHOD_BEGIN;

  GeometryLine2D geometry;
  CellGeometry* geometryLD = geometry.geometryLowerDim();
  CPPUNIT_ASSERT(!geometryLD);

  PYLITH_METHOD_END;
} // testGeomLowerDim


// End of file 
