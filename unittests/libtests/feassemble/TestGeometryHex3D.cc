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

#include "TestGeometryHex3D.hh" // Implementation of class methods

#include "pylith/feassemble/GeometryHex3D.hh"
#include "pylith/feassemble/GeometryQuad3D.hh"
#include "pylith/feassemble/GeometryPoint2D.hh"

#include "data/GeomDataHex3D.hh"

// ----------------------------------------------------------------------
CPPUNIT_TEST_SUITE_REGISTRATION( pylith::feassemble::TestGeometryHex3D );

// ----------------------------------------------------------------------
// Setup test data.
void
pylith::feassemble::TestGeometryHex3D::setUp(void)
{ // setUp
  _object = new GeometryHex3D();
  _data = new GeomDataHex3D();
} // setUp

// ----------------------------------------------------------------------
// Test constructor
void
pylith::feassemble::TestGeometryHex3D::testConstructor(void)
{ // testConstructor
  GeometryHex3D geometry;
} // testConstructor

// ----------------------------------------------------------------------
// Test geometryLowerDim().
void
pylith::feassemble::TestGeometryHex3D::testGeomLowerDim(void)
{ // testGeomLowerDim
  GeometryHex3D geometry;
  CellGeometry* geometryLD = geometry.geometryLowerDim();
  GeometryQuad3D* geometryPt = dynamic_cast<GeometryQuad3D*>(geometryLD);
  CPPUNIT_ASSERT(0 != geometryPt);
  GeometryPoint2D* geometryPt2 = dynamic_cast<GeometryPoint2D*>(geometryLD);
  CPPUNIT_ASSERT(0 == geometryPt2);
  delete geometryLD; geometryLD = 0;
} // testGeomLowerDim


// End of file 
