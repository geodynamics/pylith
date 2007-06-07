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

#include "TestGeometryQuad3D.hh" // Implementation of class methods

#include "pylith/feassemble/GeometryQuad3D.hh"
#include "pylith/feassemble/GeometryLine3D.hh"
#include "pylith/feassemble/GeometryPoint2D.hh"

#include "data/GeomDataQuad3D.hh"

// ----------------------------------------------------------------------
CPPUNIT_TEST_SUITE_REGISTRATION( pylith::feassemble::TestGeometryQuad3D );

// ----------------------------------------------------------------------
// Setup test data.
void
pylith::feassemble::TestGeometryQuad3D::setUp(void)
{ // setUp
  _object = new GeometryQuad3D();
  _data = new GeomDataQuad3D();
} // setUp

// ----------------------------------------------------------------------
// Test constructor
void
pylith::feassemble::TestGeometryQuad3D::testConstructor(void)
{ // testConstructor
  GeometryQuad3D geometry;
} // testConstructor

// ----------------------------------------------------------------------
// Test geometryLowerDim().
void
pylith::feassemble::TestGeometryQuad3D::testGeomLowerDim(void)
{ // testGeomLowerDim
  GeometryQuad3D geometry;
  CellGeometry* geometryLD = geometry.geometryLowerDim();
  GeometryLine3D* geometryPt = dynamic_cast<GeometryLine3D*>(geometryLD);
  CPPUNIT_ASSERT(0 != geometryPt);
  GeometryPoint2D* geometryPt2 = dynamic_cast<GeometryPoint2D*>(geometryLD);
  CPPUNIT_ASSERT(0 == geometryPt2);
  delete geometryLD; geometryLD = 0;
} // testGeomLowerDim


// End of file 
