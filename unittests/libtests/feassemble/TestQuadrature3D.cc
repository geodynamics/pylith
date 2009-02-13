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
#include "pylith/feassemble/QuadratureRefCell.hh"
#include "pylith/feassemble/GeometryTet3D.hh"

#include "data/QuadratureData3DLinear.hh"
#include "data/QuadratureData3DQuadratic.hh"

// ----------------------------------------------------------------------
CPPUNIT_TEST_SUITE_REGISTRATION( pylith::feassemble::TestQuadrature3D );

// ----------------------------------------------------------------------
// Test constructor
void
pylith::feassemble::TestQuadrature3D::testConstructor(void)
{ // testConstructor
  QuadratureRefCell refCell;
  Quadrature3D quadrature(refCell);
} // testConstructor

// ----------------------------------------------------------------------
// Test computeGeometry() w/linear basis fns
void
pylith::feassemble::TestQuadrature3D::testLinear(void)
{ // testLinear
  GeometryTet3D geometry;
  QuadratureRefCell refCell;
  refCell.refGeometry(&geometry);

  Quadrature3D q(refCell);
  QuadratureData3DLinear data;

  _testComputeGeometry(&q, &refCell, data);
} // testLinear

// ----------------------------------------------------------------------
// Test computeGeometry() w/quadratic basis fns
void
pylith::feassemble::TestQuadrature3D::testQuadratic(void)
{ // testQuadratic
  GeometryTet3D geometry;
  QuadratureRefCell refCell;
  refCell.refGeometry(&geometry);

  Quadrature3D q(refCell);
  QuadratureData3DQuadratic data;

  _testComputeGeometry(&q, &refCell, data);
} // testQuadratic


// End of file 
