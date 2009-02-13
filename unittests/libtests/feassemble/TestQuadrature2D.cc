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
#include "pylith/feassemble/QuadratureRefCell.hh"
#include "pylith/feassemble/GeometryTri2D.hh"

#include "data/QuadratureData2DLinear.hh"
#include "data/QuadratureData2DQuadratic.hh"

// ----------------------------------------------------------------------
CPPUNIT_TEST_SUITE_REGISTRATION( pylith::feassemble::TestQuadrature2D );

// ----------------------------------------------------------------------
// Test constructor
void
pylith::feassemble::TestQuadrature2D::testConstructor(void)
{ // testConstructor
  QuadratureRefCell refCell;
  Quadrature2D quadrature(refCell);
} // testConstructor

// ----------------------------------------------------------------------
// Test computeGeometry() w/linear basis fns
void
pylith::feassemble::TestQuadrature2D::testLinear(void)
{ // testLinear
  GeometryTri2D geometry;
  QuadratureRefCell refCell;
  refCell.refGeometry(&geometry);

  Quadrature2D q(refCell);
  QuadratureData2DLinear data;

  _testComputeGeometry(&q, &refCell, data);
} // testLinear

// ----------------------------------------------------------------------
// Test computeGeometry() w/quadratic basis fns
void
pylith::feassemble::TestQuadrature2D::testQuadratic(void)
{ // testQuadratic
  GeometryTri2D geometry;
  QuadratureRefCell refCell;
  refCell.refGeometry(&geometry);

  Quadrature2D q(refCell);
  QuadratureData2DQuadratic data;

  _testComputeGeometry(&q, &refCell, data);
} // testQuadratic


// End of file 
