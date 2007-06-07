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

#include "TestGeometryTet3D.hh" // Implementation of class methods

#include "pylith/feassemble/GeometryTet3D.hh"
#include "pylith/feassemble/GeometryTri3D.hh"
#include "pylith/feassemble/GeometryPoint2D.hh"

#include "data/GeomDataTet3D.hh"

// ----------------------------------------------------------------------
CPPUNIT_TEST_SUITE_REGISTRATION( pylith::feassemble::TestGeometryTet3D );

// ----------------------------------------------------------------------
// Setup test data.
void
pylith::feassemble::TestGeometryTet3D::setUp(void)
{ // setUp
  _object = new GeometryTet3D();
  _data = new GeomDataTet3D();
} // setUp

// ----------------------------------------------------------------------
// Test constructor
void
pylith::feassemble::TestGeometryTet3D::testConstructor(void)
{ // testConstructor
  GeometryTet3D geometry;
} // testConstructor

// ----------------------------------------------------------------------
// Test geometryLowerDim().
void
pylith::feassemble::TestGeometryTet3D::testGeomLowerDim(void)
{ // testGeomLowerDim
  GeometryTet3D geometry;
  CellGeometry* geometryLD = geometry.geometryLowerDim();
  GeometryTri3D* geometryPt = dynamic_cast<GeometryTri3D*>(geometryLD);
  CPPUNIT_ASSERT(0 != geometryPt);
  GeometryPoint2D* geometryPt2 = dynamic_cast<GeometryPoint2D*>(geometryLD);
  CPPUNIT_ASSERT(0 == geometryPt2);
  delete geometryLD; geometryLD = 0;
} // testGeomLowerDim


// End of file 
