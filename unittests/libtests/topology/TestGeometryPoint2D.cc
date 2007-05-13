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

#include "TestGeometryPoint2D.hh" // Implementation of class methods

#include "pylith/topology/GeometryPoint2D.hh"
#include "pylith/utils/array.hh" // USES double_array

#include "data/GeomDataPoint2D.hh"

// ----------------------------------------------------------------------
CPPUNIT_TEST_SUITE_REGISTRATION( pylith::topology::TestGeometryPoint2D );

// ----------------------------------------------------------------------
// Test constructor
void
pylith::topology::TestGeometryPoint2D::testConstructor(void)
{ // testConstructor
  GeometryPoint2D geometry;
} // testConstructor

// ----------------------------------------------------------------------
// Test cellDim()
void
pylith::topology::TestGeometryPoint2D::testCellDim(void)
{ // testCellDim
  GeometryPoint2D geometry;
  GeomDataPoint2D data;
  _testCellDim(geometry, data);
} // testCellDim

// ----------------------------------------------------------------------
// Test spaceDim()
void
pylith::topology::TestGeometryPoint2D::testSpaceDim(void)
{ // testSpaceDim
  GeometryPoint2D geometry;
  GeomDataPoint2D data;
  _testSpaceDim(geometry, data);
} // testSpaceDim

// ----------------------------------------------------------------------
// Test numCorners()
void
pylith::topology::TestGeometryPoint2D::testNumCorners(void)
{ // testNumCorners
  GeometryPoint2D geometry;
  GeomDataPoint2D data;
  _testNumCorners(geometry, data);
} // testNumCorners

// ----------------------------------------------------------------------
// Test geometryLowerDim().
void
pylith::topology::TestGeometryPoint2D::testGeomLowerDim(void)
{ // testGeomLowerDim
  GeometryPoint2D geometry;
  CellGeometry* geometryLD = geometry.geometryLowerDim();
  CPPUNIT_ASSERT(0 == geometryLD);
  delete geometryLD; geometryLD = 0;
} // testGeomLowerDim

// ----------------------------------------------------------------------
// Test jacobian().
void
pylith::topology::TestGeometryPoint2D::testJacobian(void)
{ // testJacobian
  GeometryPoint2D geometry;
  GeomDataPoint2D data;

  const int cellDim = data.cellDim;
  const int spaceDim = data.spaceDim;
  const int numCorners = data.numCorners;

  const int numLocs = data.numLocs;

  CPPUNIT_ASSERT_EQUAL(cellDim, geometry.cellDim());
  CPPUNIT_ASSERT_EQUAL(spaceDim, geometry.spaceDim());
  CPPUNIT_ASSERT_EQUAL(numCorners, geometry.numCorners());

  double_array jacobian(1);
  for (int iLoc=0; iLoc < numLocs; ++iLoc) {
    double_array vertices(data.vertices, numCorners*spaceDim);
    double_array location(&data.locations[iLoc], 1);

    geometry.jacobian(&jacobian, vertices, location);
    CPPUNIT_ASSERT_EQUAL(1.0, jacobian[0]);
  } //for
} // testJacobian


// End of file 
