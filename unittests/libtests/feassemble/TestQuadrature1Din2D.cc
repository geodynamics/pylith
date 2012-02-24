// -*- C++ -*-
//
// ----------------------------------------------------------------------
//
// Brad T. Aagaard, U.S. Geological Survey
// Charles A. Williams, GNS Science
// Matthew G. Knepley, University of Chicago
//
// This code was developed as part of the Computational Infrastructure
// for Geodynamics (http://geodynamics.org).
//
// Copyright (c) 2010-2012 University of California, Davis
//
// See COPYING for license information.
//
// ----------------------------------------------------------------------
//

#include <portinfo>

#include "TestQuadrature1Din2D.hh" // Implementation of class methods

#include "pylith/feassemble/Quadrature1Din2D.hh"
#include "pylith/feassemble/QuadratureRefCell.hh"
#include "pylith/feassemble/GeometryLine2D.hh"

#include "data/QuadratureData1Din2DLinear.hh"
#include "data/QuadratureData1Din2DQuadratic.hh"

// ----------------------------------------------------------------------
CPPUNIT_TEST_SUITE_REGISTRATION( pylith::feassemble::TestQuadrature1Din2D );

// ----------------------------------------------------------------------
// Test constructor
void
pylith::feassemble::TestQuadrature1Din2D::testConstructor(void)
{ // testConstructor
  QuadratureRefCell refCell;
  Quadrature1Din2D quadrature(refCell);
} // testConstructor

// ----------------------------------------------------------------------
// Test computeGeometry() w/linear basis fns
void
pylith::feassemble::TestQuadrature1Din2D::testLinear(void)
{ // testLinear
  GeometryLine2D geometry;
  QuadratureRefCell refCell;
  refCell.refGeometry(&geometry);

  Quadrature1Din2D q(refCell);
  QuadratureData1Din2DLinear data;

  _testComputeGeometry(&q, &refCell, data);
} // testLinear

// ----------------------------------------------------------------------
// Test computeGeometry() w/quadratic basis fns
void
pylith::feassemble::TestQuadrature1Din2D::testQuadratic(void)
{ // testQuadratic
  GeometryLine2D geometry;
  QuadratureRefCell refCell;
  refCell.refGeometry(&geometry);

  Quadrature1Din2D q(refCell);
  QuadratureData1Din2DQuadratic data;

  _testComputeGeometry(&q, &refCell, data);
} // testQuadratic

// End of file 
