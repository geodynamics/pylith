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
#include "pylith/feassemble/QuadratureRefCell.hh"
#include "pylith/feassemble/GeometryLine1D.hh"

#include "data/QuadratureData1DLinear.hh"
#include "data/QuadratureData1DQuadratic.hh"

// ----------------------------------------------------------------------
CPPUNIT_TEST_SUITE_REGISTRATION( pylith::feassemble::TestQuadrature1D );

// ----------------------------------------------------------------------
// Test constructor
void
pylith::feassemble::TestQuadrature1D::testConstructor(void)
{ // testConstructor
  QuadratureRefCell refCell;
  Quadrature1D quadrature(refCell);
} // testConstructor

// ----------------------------------------------------------------------
// Test computeGeometry() w/linear basis fns
void
pylith::feassemble::TestQuadrature1D::testLinear(void)
{ // testLinear
  GeometryLine1D geometry;
  QuadratureRefCell refCell;
  refCell.refGeometry(&geometry);

  Quadrature1D q(refCell);
  QuadratureData1DLinear data;

  _testComputeGeometry(&q, &refCell, data);
} // testLinear

// ----------------------------------------------------------------------
// Test computeGeometry() w/quadratic basis fns
void
pylith::feassemble::TestQuadrature1D::testQuadratic(void)
{ // testQuadratic
  GeometryLine1D geometry;
  QuadratureRefCell refCell;
  refCell.refGeometry(&geometry);

  Quadrature1D q(refCell);
  QuadratureData1DQuadratic data;

  _testComputeGeometry(&q, &refCell, data);
} // testQuadratic


// End of file 
