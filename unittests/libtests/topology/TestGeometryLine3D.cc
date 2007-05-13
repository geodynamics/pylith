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

#include "TestGeometryLine3D.hh" // Implementation of class methods

#include "pylith/topology/GeometryLine3D.hh"
#include "pylith/topology/GeometryPoint3D.hh"

#include "pylith/utils/array.hh" // USES double_array

#include "data/GeomDataLine3D.hh"

// ----------------------------------------------------------------------
CPPUNIT_TEST_SUITE_REGISTRATION( pylith::topology::TestGeometryLine3D );

// ----------------------------------------------------------------------
// Test constructor
void
pylith::topology::TestGeometryLine3D::testConstructor(void)
{ // testConstructor
  GeometryLine3D geometry;
} // testConstructor

// ----------------------------------------------------------------------
// Test cellDim()
void
pylith::topology::TestGeometryLine3D::testCellDim(void)
{ // testCellDim
  GeometryLine3D geometry;
  GeomDataLine3D data;
  _testCellDim(geometry, data);
} // testCellDim

// ----------------------------------------------------------------------
// Test spaceDim()
void
pylith::topology::TestGeometryLine3D::testSpaceDim(void)
{ // testSpaceDim
  GeometryLine3D geometry;
  GeomDataLine3D data;
  _testSpaceDim(geometry, data);
} // testSpaceDim

// ----------------------------------------------------------------------
// Test numCorners()
void
pylith::topology::TestGeometryLine3D::testNumCorners(void)
{ // testNumCorners
  GeometryLine3D geometry;
  GeomDataLine3D data;
  _testNumCorners(geometry, data);
} // testNumCorners

// ----------------------------------------------------------------------
// Test geometryLowerDim().
void
pylith::topology::TestGeometryLine3D::testGeomLowerDim(void)
{ // testGeomLowerDim
  GeometryLine3D geometry;
  CellGeometry* geometryLD = geometry.geometryLowerDim();
  GeometryPoint3D* geometryPt = dynamic_cast<GeometryPoint3D*>(geometryLD);
  CPPUNIT_ASSERT(0 != geometryLD);
  delete geometryLD; geometryLD = 0;
} // testGeomLowerDim

// ----------------------------------------------------------------------
// Test jacobian().
void
pylith::topology::TestGeometryLine3D::testJacobian(void)
{ // testJacobian
  GeometryLine3D geometry;
  GeomDataLine3D data;

  _testJacobian(&geometry, data);
} // testJacobian


// End of file 
