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

#include "TestGeometryQuad2D.hh" // Implementation of class methods

#include "pylith/feassemble/GeometryQuad2D.hh"
#include "pylith/feassemble/GeometryLine2D.hh"
#include "pylith/feassemble/GeometryPoint2D.hh"

#include "data/GeomDataQuad2D.hh"

// ----------------------------------------------------------------------
CPPUNIT_TEST_SUITE_REGISTRATION( pylith::feassemble::TestGeometryQuad2D );

// ----------------------------------------------------------------------
// Setup test data.
void
pylith::feassemble::TestGeometryQuad2D::setUp(void)
{ // setUp
  _object = new GeometryQuad2D();
  _data = new GeomDataQuad2D();
} // setUp

// ----------------------------------------------------------------------
// Test constructor
void
pylith::feassemble::TestGeometryQuad2D::testConstructor(void)
{ // testConstructor
  GeometryQuad2D geometry;
} // testConstructor

// ----------------------------------------------------------------------
// Test geometryLowerDim().
void
pylith::feassemble::TestGeometryQuad2D::testGeomLowerDim(void)
{ // testGeomLowerDim
  GeometryQuad2D geometry;
  CellGeometry* geometryLD = geometry.geometryLowerDim();
  GeometryLine2D* geometryPt = dynamic_cast<GeometryLine2D*>(geometryLD);
  CPPUNIT_ASSERT(0 != geometryPt);
  GeometryPoint2D* geometryPt2 = dynamic_cast<GeometryPoint2D*>(geometryLD);
  CPPUNIT_ASSERT(0 == geometryPt2);
  delete geometryLD; geometryLD = 0;
} // testGeomLowerDim


// End of file 
