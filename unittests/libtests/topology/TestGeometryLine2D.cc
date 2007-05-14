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

#include "pylith/topology/GeometryLine2D.hh"
#include "pylith/topology/GeometryPoint2D.hh"
#include "pylith/topology/GeometryPoint1D.hh"

#include "pylith/utils/array.hh" // USES double_array

#include "data/GeomDataLine2D.hh"

// ----------------------------------------------------------------------
CPPUNIT_TEST_SUITE_REGISTRATION( pylith::topology::TestGeometryLine2D );

// ----------------------------------------------------------------------
// Test constructor
void
pylith::topology::TestGeometryLine2D::testConstructor(void)
{ // testConstructor
  GeometryLine2D geometry;
} // testConstructor

// ----------------------------------------------------------------------
// Test cellDim()
void
pylith::topology::TestGeometryLine2D::testCellDim(void)
{ // testCellDim
  GeometryLine2D geometry;
  GeomDataLine2D data;
  _testCellDim(geometry, data);
} // testCellDim

// ----------------------------------------------------------------------
// Test spaceDim()
void
pylith::topology::TestGeometryLine2D::testSpaceDim(void)
{ // testSpaceDim
  GeometryLine2D geometry;
  GeomDataLine2D data;
  _testSpaceDim(geometry, data);
} // testSpaceDim

// ----------------------------------------------------------------------
// Test numCorners()
void
pylith::topology::TestGeometryLine2D::testNumCorners(void)
{ // testNumCorners
  GeometryLine2D geometry;
  GeomDataLine2D data;
  _testNumCorners(geometry, data);
} // testNumCorners

// ----------------------------------------------------------------------
// Test geometryLowerDim().
void
pylith::topology::TestGeometryLine2D::testGeomLowerDim(void)
{ // testGeomLowerDim
  GeometryLine2D geometry;
  CellGeometry* geometryLD = geometry.geometryLowerDim();
  GeometryPoint2D* geometryPt = dynamic_cast<GeometryPoint2D*>(geometryLD);
  CPPUNIT_ASSERT(0 != geometryPt);
  GeometryPoint1D* geometryPt2 = dynamic_cast<GeometryPoint1D*>(geometryLD);
  CPPUNIT_ASSERT(0 == geometryPt2);
  delete geometryLD; geometryLD = 0;
} // testGeomLowerDim

// ----------------------------------------------------------------------
// Test jacobian().
void
pylith::topology::TestGeometryLine2D::testJacobian(void)
{ // testJacobian
  GeometryLine2D geometry;
  GeomDataLine2D data;

  _testJacobian(&geometry, data);
} // testJacobian


// End of file 
