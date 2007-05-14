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

#include "pylith/topology/GeometryQuad2D.hh"
#include "pylith/topology/GeometryLine2D.hh"
#include "pylith/topology/GeometryPoint2D.hh"

#include "pylith/utils/array.hh" // USES double_array

#include "data/GeomDataQuad2D.hh"

// ----------------------------------------------------------------------
CPPUNIT_TEST_SUITE_REGISTRATION( pylith::topology::TestGeometryQuad2D );

// ----------------------------------------------------------------------
// Test constructor
void
pylith::topology::TestGeometryQuad2D::testConstructor(void)
{ // testConstructor
  GeometryQuad2D geometry;
} // testConstructor

// ----------------------------------------------------------------------
// Test cellDim()
void
pylith::topology::TestGeometryQuad2D::testCellDim(void)
{ // testCellDim
  GeometryQuad2D geometry;
  GeomDataQuad2D data;
  _testCellDim(geometry, data);
} // testCellDim

// ----------------------------------------------------------------------
// Test spaceDim()
void
pylith::topology::TestGeometryQuad2D::testSpaceDim(void)
{ // testSpaceDim
  GeometryQuad2D geometry;
  GeomDataQuad2D data;
  _testSpaceDim(geometry, data);
} // testSpaceDim

// ----------------------------------------------------------------------
// Test numCorners()
void
pylith::topology::TestGeometryQuad2D::testNumCorners(void)
{ // testNumCorners
  GeometryQuad2D geometry;
  GeomDataQuad2D data;
  _testNumCorners(geometry, data);
} // testNumCorners

// ----------------------------------------------------------------------
// Test geometryLowerDim().
void
pylith::topology::TestGeometryQuad2D::testGeomLowerDim(void)
{ // testGeomLowerDim
  GeometryQuad2D geometry;
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
pylith::topology::TestGeometryQuad2D::testJacobian(void)
{ // testJacobian
  GeometryQuad2D geometry;
  GeomDataQuad2D data;

  _testJacobian(&geometry, data);
} // testJacobian


// End of file 
