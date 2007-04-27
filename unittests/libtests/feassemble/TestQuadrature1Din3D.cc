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

#include "TestQuadrature1Din3D.hh" // Implementation of class methods

#include "pylith/feassemble/Quadrature1Din3D.hh"

#include "data/QuadratureData1Din3DLinear.hh"
#include "data/QuadratureData1Din3DQuadratic.hh"

// ----------------------------------------------------------------------
CPPUNIT_TEST_SUITE_REGISTRATION( pylith::feassemble::TestQuadrature1Din3D );

// ----------------------------------------------------------------------
// Test constructor
void
pylith::feassemble::TestQuadrature1Din3D::testConstructor(void)
{ // testConstructor
  Quadrature1Din3D quadrature;
} // testConstructor

// ----------------------------------------------------------------------
// Test computeGeometry() w/linear basis fns
void
pylith::feassemble::TestQuadrature1Din3D::testLinear(void)
{ // testLinear
  Quadrature1Din3D q;
  QuadratureData1Din3DLinear data;

  _testComputeGeometryVert(&q, data);
  _testComputeGeometryQuad(&q, data);
} // testLinear

// ----------------------------------------------------------------------
// Test computeGeometry() w/quadratic basis fns
void
pylith::feassemble::TestQuadrature1Din3D::testQuadratic(void)
{ // testQuadratic
  Quadrature1Din3D q;
  QuadratureData1Din3DQuadratic data;

  _testComputeGeometryVert(&q, data);
  _testComputeGeometryQuad(&q, data);
} // testQuadratic

// End of file 
