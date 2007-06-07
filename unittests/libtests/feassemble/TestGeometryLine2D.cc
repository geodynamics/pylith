// -*- C++ -*-
//
// ----------------------------------------------------------------------
//
//                           Brad T. Aagaard
//                        U.S. Geological Survey
//
// {LicenseText}
//
// ----------------------------------------------------------------------
//

#include <portinfo>

#include "TestGeometryLine2D.hh" // Implementation of class methods

#include "pylith/feassemble/GeometryLine2D.hh"
#include "pylith/feassemble/GeometryPoint2D.hh"
#include "pylith/feassemble/GeometryPoint1D.hh"

#include "data/GeomDataLine2D.hh"

// ----------------------------------------------------------------------
CPPUNIT_TEST_SUITE_REGISTRATION( pylith::feassemble::TestGeometryLine2D );

// ----------------------------------------------------------------------
// Setup test data.
void
pylith::feassemble::TestGeometryLine2D::setUp(void)
{ // setUp
  _object = new GeometryLine2D();
  _data = new GeomDataLine2D();
} // setUp

// ----------------------------------------------------------------------
// Test constructor
void
pylith::feassemble::TestGeometryLine2D::testConstructor(void)
{ // testConstructor
  GeometryLine2D geometry;
} // testConstructor

// ----------------------------------------------------------------------
// Test geometryLowerDim().
void
pylith::feassemble::TestGeometryLine2D::testGeomLowerDim(void)
{ // testGeomLowerDim
  GeometryLine2D geometry;
  CellGeometry* geometryLD = geometry.geometryLowerDim();
  GeometryPoint2D* geometryPt = dynamic_cast<GeometryPoint2D*>(geometryLD);
  CPPUNIT_ASSERT(0 != geometryPt);
  GeometryPoint1D* geometryPt2 = dynamic_cast<GeometryPoint1D*>(geometryLD);
  CPPUNIT_ASSERT(0 == geometryPt2);
  delete geometryLD; geometryLD = 0;
} // testGeomLowerDim


// End of file 
