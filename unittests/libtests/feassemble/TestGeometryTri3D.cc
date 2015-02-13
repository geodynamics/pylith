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

#include "TestGeometryTri3D.hh" // Implementation of class methods

#include "pylith/feassemble/GeometryTri3D.hh"
#include "pylith/feassemble/GeometryLine3D.hh"
#include "pylith/feassemble/GeometryLine2D.hh"

#include "data/GeomDataTri3D.hh"

#include "pylith/utils/error.h" // USES PYLITH_METHOD_BEGIN/END

// ----------------------------------------------------------------------
CPPUNIT_TEST_SUITE_REGISTRATION( pylith::feassemble::TestGeometryTri3D );

// ----------------------------------------------------------------------
// Setup test data.
void
pylith::feassemble::TestGeometryTri3D::setUp(void)
{ // setUp
  PYLITH_METHOD_BEGIN;

  _object = new GeometryTri3D();
  _data = new GeomDataTri3D();

  PYLITH_METHOD_END;
} // setUp

// ----------------------------------------------------------------------
// Test constructor
void
pylith::feassemble::TestGeometryTri3D::testConstructor(void)
{ // testConstructor
  PYLITH_METHOD_BEGIN;

  GeometryTri3D geometry;

  PYLITH_METHOD_END;
} // testConstructor

// ----------------------------------------------------------------------
// Test geometryLowerDim().
void
pylith::feassemble::TestGeometryTri3D::testGeomLowerDim(void)
{ // testGeomLowerDim
  PYLITH_METHOD_BEGIN;

  GeometryTri3D geometry;
  CellGeometry* geometryLD = geometry.geometryLowerDim();
  GeometryLine3D* geometryPt = dynamic_cast<GeometryLine3D*>(geometryLD);
  CPPUNIT_ASSERT(geometryPt);
  GeometryLine2D* geometryPt2 = dynamic_cast<GeometryLine2D*>(geometryLD);
  CPPUNIT_ASSERT(!geometryPt2);
  delete geometryLD; geometryLD = 0;

  PYLITH_METHOD_END;
} // testGeomLowerDim


// End of file 
