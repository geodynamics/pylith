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

#include "TestQuadrature1D.hh" // Implementation of class methods

#include "pylith/feassemble/Quadrature1D.hh"

#include "data/QuadratureData1DLinear.hh"
#include "data/QuadratureData1DQuadratic.hh"

// ----------------------------------------------------------------------
CPPUNIT_TEST_SUITE_REGISTRATION( pylith::feassemble::TestQuadrature1D );

// ----------------------------------------------------------------------
// Test constructor
void
pylith::feassemble::TestQuadrature1D::testConstructor(void)
{ // testConstructor
  Quadrature1D quadrature;
} // testConstructor

// ----------------------------------------------------------------------
// Test computeGeometry() w/linear basis fns
void
pylith::feassemble::TestQuadrature1D::testLinear(void)
{ // testLinear
  Quadrature1D q;
  QuadratureData1DLinear data;

  _testComputeGeometryVert(&q, data);
  _testComputeGeometryQuad(&q, data);
} // testLinear

// ----------------------------------------------------------------------
// Test computeGeometry() w/quadratic basis fns
void
pylith::feassemble::TestQuadrature1D::testQuadratic(void)
{ // testQuadratic
  Quadrature1D q;
  QuadratureData1DQuadratic data;

  _testComputeGeometryVert(&q, data);
  _testComputeGeometryQuad(&q, data);
} // testQuadratic


// End of file 
