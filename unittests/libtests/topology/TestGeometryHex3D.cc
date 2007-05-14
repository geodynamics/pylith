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

#include "pylith/topology/GeometryHex3D.hh"
#include "pylith/topology/GeometryQuad3D.hh"
#include "pylith/topology/GeometryPoint2D.hh"

#include "pylith/utils/array.hh" // USES double_array

#include "data/GeomDataHex3D.hh"

// ----------------------------------------------------------------------
CPPUNIT_TEST_SUITE_REGISTRATION( pylith::topology::TestGeometryHex3D );

// ----------------------------------------------------------------------
// Test constructor
void
pylith::topology::TestGeometryHex3D::testConstructor(void)
{ // testConstructor
  GeometryHex3D geometry;
} // testConstructor

// ----------------------------------------------------------------------
// Test cellDim()
void
pylith::topology::TestGeometryHex3D::testCellDim(void)
{ // testCellDim
  GeometryHex3D geometry;
  GeomDataHex3D data;
  _testCellDim(geometry, data);
} // testCellDim

// ----------------------------------------------------------------------
// Test spaceDim()
void
pylith::topology::TestGeometryHex3D::testSpaceDim(void)
{ // testSpaceDim
  GeometryHex3D geometry;
  GeomDataHex3D data;
  _testSpaceDim(geometry, data);
} // testSpaceDim

// ----------------------------------------------------------------------
// Test numCorners()
void
pylith::topology::TestGeometryHex3D::testNumCorners(void)
{ // testNumCorners
  GeometryHex3D geometry;
  GeomDataHex3D data;
  _testNumCorners(geometry, data);
} // testNumCorners

// ----------------------------------------------------------------------
// Test geometryLowerDim().
void
pylith::topology::TestGeometryHex3D::testGeomLowerDim(void)
{ // testGeomLowerDim
  GeometryHex3D geometry;
  CellGeometry* geometryLD = geometry.geometryLowerDim();
  GeometryQuad3D* geometryPt = dynamic_cast<GeometryQuad3D*>(geometryLD);
  CPPUNIT_ASSERT(0 != geometryPt);
  GeometryPoint2D* geometryPt2 = dynamic_cast<GeometryPoint2D*>(geometryLD);
  CPPUNIT_ASSERT(0 == geometryPt2);
  delete geometryLD; geometryLD = 0;
} // testGeomLowerDim

// ----------------------------------------------------------------------
// Test jacobian().
void
pylith::topology::TestGeometryHex3D::testJacobian(void)
{ // testJacobian
  GeometryHex3D geometry;
  GeomDataHex3D data;

  _testJacobian(&geometry, data);
} // testJacobian


// End of file 
