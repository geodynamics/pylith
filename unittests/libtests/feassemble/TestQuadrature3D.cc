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

#include "TestQuadrature3D.hh" // Implementation of class methods

#include "pylith/feassemble/Quadrature3D.hh"

#include "data/QuadratureData3DLinear.hh"
#include "data/QuadratureData3DQuadratic.hh"

// ----------------------------------------------------------------------
CPPUNIT_TEST_SUITE_REGISTRATION( pylith::feassemble::TestQuadrature3D );

// ----------------------------------------------------------------------
// Test constructor
void
pylith::feassemble::TestQuadrature3D::testConstructor(void)
{ // testConstructor
  Quadrature3D quadrature;
} // testConstructor

// ----------------------------------------------------------------------
// Test computeGeometry() w/linear basis fns
void
pylith::feassemble::TestQuadrature3D::testLinear(void)
{ // testLinear
  Quadrature3D q;
  QuadratureData3DLinear data;

  _testComputeGeometry(&q, data);
} // testLinear

// ----------------------------------------------------------------------
// Test computeGeometry() w/quadratic basis fns
void
pylith::feassemble::TestQuadrature3D::testQuadratic(void)
{ // testQuadratic
  Quadrature3D q;
  QuadratureData3DQuadratic data;

  _testComputeGeometry(&q, data);
} // testQuadratic

// End of file 
