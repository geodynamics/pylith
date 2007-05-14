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

#include "TestGeometryTri2D.hh" // Implementation of class methods

#include "pylith/topology/GeometryTri2D.hh"
#include "pylith/topology/GeometryLine2D.hh"
#include "pylith/topology/GeometryPoint2D.hh"

#include "pylith/utils/array.hh" // USES double_array

#include "data/GeomDataTri2D.hh"

// ----------------------------------------------------------------------
CPPUNIT_TEST_SUITE_REGISTRATION( pylith::topology::TestGeometryTri2D );

// ----------------------------------------------------------------------
// Test constructor
void
pylith::topology::TestGeometryTri2D::testConstructor(void)
{ // testConstructor
  GeometryTri2D geometry;
} // testConstructor

// ----------------------------------------------------------------------
// Test cellDim()
void
pylith::topology::TestGeometryTri2D::testCellDim(void)
{ // testCellDim
  GeometryTri2D geometry;
  GeomDataTri2D data;
  _testCellDim(geometry, data);
} // testCellDim

// ----------------------------------------------------------------------
// Test spaceDim()
void
pylith::topology::TestGeometryTri2D::testSpaceDim(void)
{ // testSpaceDim
  GeometryTri2D geometry;
  GeomDataTri2D data;
  _testSpaceDim(geometry, data);
} // testSpaceDim

// ----------------------------------------------------------------------
// Test numCorners()
void
pylith::topology::TestGeometryTri2D::testNumCorners(void)
{ // testNumCorners
  GeometryTri2D geometry;
  GeomDataTri2D data;
  _testNumCorners(geometry, data);
} // testNumCorners

// ----------------------------------------------------------------------
// Test geometryLowerDim().
void
pylith::topology::TestGeometryTri2D::testGeomLowerDim(void)
{ // testGeomLowerDim
  GeometryTri2D geometry;
  CellGeometry* geometryLD = geometry.geometryLowerDim();
  GeometryLine2D* geometryPt = dynamic_cast<GeometryLine2D*>(geometryLD);
  CPPUNIT_ASSERT(0 != geometryPt);
  GeometryPoint2D* geometryPt2 = dynamic_cast<GeometryPoint2D*>(geometryLD);
  CPPUNIT_ASSERT(0 == geometryPt2);
  delete geometryLD; geometryLD = 0;
} // testGeomLowerDim

// ----------------------------------------------------------------------
// Test jacobian().
void
pylith::topology::TestGeometryTri2D::testJacobian(void)
{ // testJacobian
  GeometryTri2D geometry;
  GeomDataTri2D data;

  _testJacobian(&geometry, data);
} // testJacobian


// End of file 
