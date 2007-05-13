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

#include "TestGeometryPoint1D.hh" // Implementation of class methods

#include "pylith/topology/GeometryPoint1D.hh"
#include "pylith/utils/array.hh" // USES double_array

#include "data/GeomDataPoint1D.hh"

// ----------------------------------------------------------------------
CPPUNIT_TEST_SUITE_REGISTRATION( pylith::topology::TestGeometryPoint1D );

// ----------------------------------------------------------------------
// Test constructor
void
pylith::topology::TestGeometryPoint1D::testConstructor(void)
{ // testConstructor
  GeometryPoint1D geometry;
} // testConstructor

// ----------------------------------------------------------------------
// Test cellDim()
void
pylith::topology::TestGeometryPoint1D::testCellDim(void)
{ // testCellDim
  GeometryPoint1D geometry;
  GeomDataPoint1D data;
  _testCellDim(geometry, data);
} // testCellDim

// ----------------------------------------------------------------------
// Test spaceDim()
void
pylith::topology::TestGeometryPoint1D::testSpaceDim(void)
{ // testSpaceDim
  GeometryPoint1D geometry;
  GeomDataPoint1D data;
  _testSpaceDim(geometry, data);
} // testSpaceDim

// ----------------------------------------------------------------------
// Test numCorners()
void
pylith::topology::TestGeometryPoint1D::testNumCorners(void)
{ // testNumCorners
  GeometryPoint1D geometry;
  GeomDataPoint1D data;
  _testNumCorners(geometry, data);
} // testNumCorners

// ----------------------------------------------------------------------
// Test geometryLowerDim().
void
pylith::topology::TestGeometryPoint1D::testGeomLowerDim(void)
{ // testGeomLowerDim
  GeometryPoint1D geometry;
  CellGeometry* geometryLD = geometry.geometryLowerDim();
  CPPUNIT_ASSERT(0 == geometryLD);
  delete geometryLD; geometryLD = 0;
} // testGeomLowerDim

// ----------------------------------------------------------------------
// Test jacobian().
void
pylith::topology::TestGeometryPoint1D::testJacobian(void)
{ // testJacobian
  GeometryPoint1D geometry;
  GeomDataPoint1D data;

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
