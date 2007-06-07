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

#include "TestGeometryTri3D.hh" // Implementation of class methods

#include "pylith/feassemble/GeometryTri3D.hh"
#include "pylith/feassemble/GeometryLine3D.hh"
#include "pylith/feassemble/GeometryPoint2D.hh"

#include "data/GeomDataTri3D.hh"

// ----------------------------------------------------------------------
CPPUNIT_TEST_SUITE_REGISTRATION( pylith::feassemble::TestGeometryTri3D );

// ----------------------------------------------------------------------
// Setup test data.
void
pylith::feassemble::TestGeometryTri3D::setUp(void)
{ // setUp
  _object = new GeometryTri3D();
  _data = new GeomDataTri3D();
} // setUp

// ----------------------------------------------------------------------
// Test constructor
void
pylith::feassemble::TestGeometryTri3D::testConstructor(void)
{ // testConstructor
  GeometryTri3D geometry;
} // testConstructor

// ----------------------------------------------------------------------
// Test geometryLowerDim().
void
pylith::feassemble::TestGeometryTri3D::testGeomLowerDim(void)
{ // testGeomLowerDim
  GeometryTri3D geometry;
  CellGeometry* geometryLD = geometry.geometryLowerDim();
  GeometryLine3D* geometryPt = dynamic_cast<GeometryLine3D*>(geometryLD);
  CPPUNIT_ASSERT(0 != geometryPt);
  GeometryPoint2D* geometryPt2 = dynamic_cast<GeometryPoint2D*>(geometryLD);
  CPPUNIT_ASSERT(0 == geometryPt2);
  delete geometryLD; geometryLD = 0;
} // testGeomLowerDim


// End of file 
