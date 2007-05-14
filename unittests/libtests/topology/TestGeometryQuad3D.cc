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

#include "pylith/topology/GeometryQuad3D.hh"
#include "pylith/topology/GeometryLine3D.hh"
#include "pylith/topology/GeometryPoint2D.hh"

#include "pylith/utils/array.hh" // USES double_array

#include "data/GeomDataQuad3D.hh"

// ----------------------------------------------------------------------
CPPUNIT_TEST_SUITE_REGISTRATION( pylith::topology::TestGeometryQuad3D );

// ----------------------------------------------------------------------
// Test constructor
void
pylith::topology::TestGeometryQuad3D::testConstructor(void)
{ // testConstructor
  GeometryQuad3D geometry;
} // testConstructor

// ----------------------------------------------------------------------
// Test cellDim()
void
pylith::topology::TestGeometryQuad3D::testCellDim(void)
{ // testCellDim
  GeometryQuad3D geometry;
  GeomDataQuad3D data;
  _testCellDim(geometry, data);
} // testCellDim

// ----------------------------------------------------------------------
// Test spaceDim()
void
pylith::topology::TestGeometryQuad3D::testSpaceDim(void)
{ // testSpaceDim
  GeometryQuad3D geometry;
  GeomDataQuad3D data;
  _testSpaceDim(geometry, data);
} // testSpaceDim

// ----------------------------------------------------------------------
// Test numCorners()
void
pylith::topology::TestGeometryQuad3D::testNumCorners(void)
{ // testNumCorners
  GeometryQuad3D geometry;
  GeomDataQuad3D data;
  _testNumCorners(geometry, data);
} // testNumCorners

// ----------------------------------------------------------------------
// Test geometryLowerDim().
void
pylith::topology::TestGeometryQuad3D::testGeomLowerDim(void)
{ // testGeomLowerDim
  GeometryQuad3D geometry;
  CellGeometry* geometryLD = geometry.geometryLowerDim();
  GeometryLine3D* geometryPt = dynamic_cast<GeometryLine3D*>(geometryLD);
  CPPUNIT_ASSERT(0 != geometryPt);
  GeometryPoint2D* geometryPt2 = dynamic_cast<GeometryPoint2D*>(geometryLD);
  CPPUNIT_ASSERT(0 == geometryPt2);
  delete geometryLD; geometryLD = 0;
} // testGeomLowerDim

// ----------------------------------------------------------------------
// Test jacobian().
void
pylith::topology::TestGeometryQuad3D::testJacobian(void)
{ // testJacobian
  GeometryQuad3D geometry;
  GeomDataQuad3D data;

  _testJacobian(&geometry, data);
} // testJacobian


// End of file 
