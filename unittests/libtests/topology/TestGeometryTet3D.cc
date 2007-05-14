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

#include "pylith/topology/GeometryTet3D.hh"
#include "pylith/topology/GeometryTri3D.hh"
#include "pylith/topology/GeometryPoint2D.hh"

#include "pylith/utils/array.hh" // USES double_array

#include "data/GeomDataTet3D.hh"

// ----------------------------------------------------------------------
CPPUNIT_TEST_SUITE_REGISTRATION( pylith::topology::TestGeometryTet3D );

// ----------------------------------------------------------------------
// Test constructor
void
pylith::topology::TestGeometryTet3D::testConstructor(void)
{ // testConstructor
  GeometryTet3D geometry;
} // testConstructor

// ----------------------------------------------------------------------
// Test cellDim()
void
pylith::topology::TestGeometryTet3D::testCellDim(void)
{ // testCellDim
  GeometryTet3D geometry;
  GeomDataTet3D data;
  _testCellDim(geometry, data);
} // testCellDim

// ----------------------------------------------------------------------
// Test spaceDim()
void
pylith::topology::TestGeometryTet3D::testSpaceDim(void)
{ // testSpaceDim
  GeometryTet3D geometry;
  GeomDataTet3D data;
  _testSpaceDim(geometry, data);
} // testSpaceDim

// ----------------------------------------------------------------------
// Test numCorners()
void
pylith::topology::TestGeometryTet3D::testNumCorners(void)
{ // testNumCorners
  GeometryTet3D geometry;
  GeomDataTet3D data;
  _testNumCorners(geometry, data);
} // testNumCorners

// ----------------------------------------------------------------------
// Test geometryLowerDim().
void
pylith::topology::TestGeometryTet3D::testGeomLowerDim(void)
{ // testGeomLowerDim
  GeometryTet3D geometry;
  CellGeometry* geometryLD = geometry.geometryLowerDim();
  GeometryTri3D* geometryPt = dynamic_cast<GeometryTri3D*>(geometryLD);
  CPPUNIT_ASSERT(0 != geometryPt);
  GeometryPoint2D* geometryPt2 = dynamic_cast<GeometryPoint2D*>(geometryLD);
  CPPUNIT_ASSERT(0 == geometryPt2);
  delete geometryLD; geometryLD = 0;
} // testGeomLowerDim

// ----------------------------------------------------------------------
// Test jacobian().
void
pylith::topology::TestGeometryTet3D::testJacobian(void)
{ // testJacobian
  GeometryTet3D geometry;
  GeomDataTet3D data;

  _testJacobian(&geometry, data);
} // testJacobian


// End of file 
