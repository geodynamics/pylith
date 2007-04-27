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

#include "TestQuadrature1Din2D.hh" // Implementation of class methods

#include "pylith/feassemble/Quadrature1Din2D.hh"

#include "data/QuadratureData1Din2DLinear.hh"
#include "data/QuadratureData1Din2DQuadratic.hh"

// ----------------------------------------------------------------------
CPPUNIT_TEST_SUITE_REGISTRATION( pylith::feassemble::TestQuadrature1Din2D );

// ----------------------------------------------------------------------
// Test constructor
void
pylith::feassemble::TestQuadrature1Din2D::testConstructor(void)
{ // testConstructor
  Quadrature1Din2D quadrature;
} // testConstructor

// ----------------------------------------------------------------------
// Test computeGeometry() w/linear basis fns
void
pylith::feassemble::TestQuadrature1Din2D::testLinear(void)
{ // testLinear
  Quadrature1Din2D q;
  QuadratureData1Din2DLinear data;

  _testComputeGeometryVert(&q, data);
  _testComputeGeometryQuad(&q, data);
} // testLinear

// ----------------------------------------------------------------------
// Test computeGeometry() w/quadratic basis fns
void
pylith::feassemble::TestQuadrature1Din2D::testQuadratic(void)
{ // testQuadratic
  Quadrature1Din2D q;
  QuadratureData1Din2DQuadratic data;

  _testComputeGeometryVert(&q, data);
  _testComputeGeometryQuad(&q, data);
} // testQuadratic

// End of file 
