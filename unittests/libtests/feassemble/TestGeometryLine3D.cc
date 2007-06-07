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

#include "TestGeometryLine3D.hh" // Implementation of class methods

#include "pylith/feassemble/GeometryLine3D.hh"
#include "pylith/feassemble/GeometryPoint3D.hh"
#include "pylith/feassemble/GeometryPoint2D.hh"

#include "data/GeomDataLine3D.hh"

// ----------------------------------------------------------------------
CPPUNIT_TEST_SUITE_REGISTRATION( pylith::feassemble::TestGeometryLine3D );

// ----------------------------------------------------------------------
// Setup test data.
void
pylith::feassemble::TestGeometryLine3D::setUp(void)
{ // setUp
  _object = new GeometryLine3D();
  _data = new GeomDataLine3D();
} // setUp

// ----------------------------------------------------------------------
// Test constructor
void
pylith::feassemble::TestGeometryLine3D::testConstructor(void)
{ // testConstructor
  GeometryLine3D geometry;
} // testConstructor

// ----------------------------------------------------------------------
// Test geometryLowerDim().
void
pylith::feassemble::TestGeometryLine3D::testGeomLowerDim(void)
{ // testGeomLowerDim
  GeometryLine3D geometry;
  CellGeometry* geometryLD = geometry.geometryLowerDim();
  GeometryPoint3D* geometryPt = dynamic_cast<GeometryPoint3D*>(geometryLD);
  CPPUNIT_ASSERT(0 != geometryPt);
  GeometryPoint2D* geometryPt2 = dynamic_cast<GeometryPoint2D*>(geometryLD);
  CPPUNIT_ASSERT(0 == geometryPt2);
  delete geometryLD; geometryLD = 0;
} // testGeomLowerDim


// End of file 
