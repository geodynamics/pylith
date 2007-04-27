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

#include "TestQuadrature2D.hh" // Implementation of class methods

#include "pylith/feassemble/Quadrature2D.hh"

#include "data/QuadratureData2DLinear.hh"
#include "data/QuadratureData2DQuadratic.hh"

// ----------------------------------------------------------------------
CPPUNIT_TEST_SUITE_REGISTRATION( pylith::feassemble::TestQuadrature2D );

// ----------------------------------------------------------------------
// Test constructor
void
pylith::feassemble::TestQuadrature2D::testConstructor(void)
{ // testConstructor
  Quadrature2D quadrature;
} // testConstructor

// ----------------------------------------------------------------------
// Test computeGeometry() w/linear basis fns
void
pylith::feassemble::TestQuadrature2D::testLinear(void)
{ // testLinear
  Quadrature2D q;
  QuadratureData2DLinear data;

  _testComputeGeometryVert(&q, data);
  _testComputeGeometryQuad(&q, data);
} // testLinear

// ----------------------------------------------------------------------
// Test computeGeometry() w/quadratic basis fns
void
pylith::feassemble::TestQuadrature2D::testQuadratic(void)
{ // testQuadratic
  Quadrature2D q;
  QuadratureData2DQuadratic data;

  _testComputeGeometryVert(&q, data);
  _testComputeGeometryQuad(&q, data);
} // testQuadratic

// End of file 
