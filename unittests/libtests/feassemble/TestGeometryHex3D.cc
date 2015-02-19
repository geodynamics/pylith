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

#include "TestGeometryHex3D.hh" // Implementation of class methods

#include "pylith/feassemble/GeometryHex3D.hh"
#include "pylith/feassemble/GeometryQuad3D.hh"
#include "pylith/feassemble/GeometryQuad2D.hh"

#include "data/GeomDataHex3D.hh"

#include "pylith/utils/error.h" // USES PYLITH_METHOD_BEGIN/END

// ----------------------------------------------------------------------
CPPUNIT_TEST_SUITE_REGISTRATION( pylith::feassemble::TestGeometryHex3D );

// ----------------------------------------------------------------------
// Setup test data.
void
pylith::feassemble::TestGeometryHex3D::setUp(void)
{ // setUp
  PYLITH_METHOD_BEGIN;

  _object = new GeometryHex3D();
  _data = new GeomDataHex3D();

  PYLITH_METHOD_END;
} // setUp

// ----------------------------------------------------------------------
// Test constructor
void
pylith::feassemble::TestGeometryHex3D::testConstructor(void)
{ // testConstructor
  PYLITH_METHOD_BEGIN;

  GeometryHex3D geometry;

  PYLITH_METHOD_END;
} // testConstructor

// ----------------------------------------------------------------------
// Test geometryLowerDim().
void
pylith::feassemble::TestGeometryHex3D::testGeomLowerDim(void)
{ // testGeomLowerDim
  PYLITH_METHOD_BEGIN;

  GeometryHex3D geometry;
  CellGeometry* geometryLD = geometry.geometryLowerDim();
  GeometryQuad3D* geometryPt = dynamic_cast<GeometryQuad3D*>(geometryLD);
  CPPUNIT_ASSERT(geometryPt);
  GeometryQuad2D* geometryPt2 = dynamic_cast<GeometryQuad2D*>(geometryLD);
  CPPUNIT_ASSERT(!geometryPt2);
  delete geometryLD; geometryLD = 0;

  PYLITH_METHOD_END;
} // testGeomLowerDim


// End of file 
