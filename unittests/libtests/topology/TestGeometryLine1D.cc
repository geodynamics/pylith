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

#include "TestGeometryLine1D.hh" // Implementation of class methods

#include "pylith/topology/GeometryLine1D.hh"
#include "pylith/topology/GeometryPoint1D.hh"

#include "pylith/utils/array.hh" // USES double_array

#include "data/GeomDataLine1D.hh"

// ----------------------------------------------------------------------
CPPUNIT_TEST_SUITE_REGISTRATION( pylith::topology::TestGeometryLine1D );

// ----------------------------------------------------------------------
// Test constructor
void
pylith::topology::TestGeometryLine1D::testConstructor(void)
{ // testConstructor
  GeometryLine1D geometry;
} // testConstructor

// ----------------------------------------------------------------------
// Test cellDim()
void
pylith::topology::TestGeometryLine1D::testCellDim(void)
{ // testCellDim
  GeometryLine1D geometry;
  GeomDataLine1D data;
  _testCellDim(geometry, data);
} // testCellDim

// ----------------------------------------------------------------------
// Test spaceDim()
void
pylith::topology::TestGeometryLine1D::testSpaceDim(void)
{ // testSpaceDim
  GeometryLine1D geometry;
  GeomDataLine1D data;
  _testSpaceDim(geometry, data);
} // testSpaceDim

// ----------------------------------------------------------------------
// Test numCorners()
void
pylith::topology::TestGeometryLine1D::testNumCorners(void)
{ // testNumCorners
  GeometryLine1D geometry;
  GeomDataLine1D data;
  _testNumCorners(geometry, data);
} // testNumCorners

// ----------------------------------------------------------------------
// Test geometryLowerDim().
void
pylith::topology::TestGeometryLine1D::testGeomLowerDim(void)
{ // testGeomLowerDim
  GeometryLine1D geometry;
  CellGeometry* geometryLD = geometry.geometryLowerDim();
  GeometryPoint1D* geometryPt = dynamic_cast<GeometryPoint1D*>(geometryLD);
  CPPUNIT_ASSERT(0 != geometryLD);
  delete geometryLD; geometryLD = 0;
} // testGeomLowerDim

// ----------------------------------------------------------------------
// Test jacobian().
void
pylith::topology::TestGeometryLine1D::testJacobian(void)
{ // testJacobian
  GeometryLine1D geometry;
  GeomDataLine1D data;

  _testJacobian(&geometry, data);
} // testJacobian


// End of file 
