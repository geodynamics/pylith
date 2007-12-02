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

#include "TestQuadrature2Din3D.hh" // Implementation of class methods

#include "pylith/feassemble/Quadrature2Din3D.hh"
#include "pylith/feassemble/GeometryTri3D.hh"

#include "data/QuadratureData2Din3DLinearXYZ.hh"
#include "data/QuadratureData2Din3DLinearXY.hh"
#include "data/QuadratureData2Din3DLinearYZ.hh"
#include "data/QuadratureData2Din3DLinearXZ.hh"
#include "data/QuadratureData2Din3DQuadratic.hh"

// ----------------------------------------------------------------------
CPPUNIT_TEST_SUITE_REGISTRATION( pylith::feassemble::TestQuadrature2Din3D );

// ----------------------------------------------------------------------
// Test constructor
void
pylith::feassemble::TestQuadrature2Din3D::testConstructor(void)
{ // testConstructor
  Quadrature2Din3D quadrature;
} // testConstructor

// ----------------------------------------------------------------------
// Test computeGeometry() w/linear basis fns
void
pylith::feassemble::TestQuadrature2Din3D::testLinearXYZ(void)
{ // testLinearXYZ
  Quadrature2Din3D q;
  GeometryTri3D geometry;
  q.refGeometry(&geometry);
  QuadratureData2Din3DLinearXYZ data;

  _testComputeGeometry(&q, data);
} // testLinearXYZ

// ----------------------------------------------------------------------
// Test computeGeometry() w/linear basis fns
void
pylith::feassemble::TestQuadrature2Din3D::testLinearXY(void)
{ // testLinearXY
  Quadrature2Din3D q;
  QuadratureData2Din3DLinearXY data;
  GeometryTri3D geometry;
  q.refGeometry(&geometry);

  _testComputeGeometry(&q, data);
} // testLinearXY

// ----------------------------------------------------------------------
// Test computeGeometry() w/linear basis fns
void
pylith::feassemble::TestQuadrature2Din3D::testLinearYZ(void)
{ // testLinearYZ
  Quadrature2Din3D q;
  GeometryTri3D geometry;
  q.refGeometry(&geometry);
  QuadratureData2Din3DLinearYZ data;

  _testComputeGeometry(&q, data);
} // testLinearYZ

// ----------------------------------------------------------------------
// Test computeGeometry() w/linear basis fns
void
pylith::feassemble::TestQuadrature2Din3D::testLinearXZ(void)
{ // testLinearXZ
  Quadrature2Din3D q;
  GeometryTri3D geometry;
  q.refGeometry(&geometry);
  QuadratureData2Din3DLinearXZ data;

  _testComputeGeometry(&q, data);
} // testLinearXZ

// ----------------------------------------------------------------------
// Test computeGeometry() w/quadratic basis fns
void
pylith::feassemble::TestQuadrature2Din3D::testQuadratic(void)
{ // testQuadratic
  Quadrature2Din3D q;
  GeometryTri3D geometry;
  q.refGeometry(&geometry);
  QuadratureData2Din3DQuadratic data;

  _testComputeGeometry(&q, data);
} // testQuadratic

// End of file 
