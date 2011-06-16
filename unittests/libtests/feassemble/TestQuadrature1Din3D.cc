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
// Copyright (c) 2010-2011 University of California, Davis
//
// See COPYING for license information.
//
// ----------------------------------------------------------------------
//

#include <portinfo>

#include "TestQuadrature1Din3D.hh" // Implementation of class methods

#include "pylith/feassemble/Quadrature1Din3D.hh"
#include "pylith/feassemble/QuadratureRefCell.hh"
#include "pylith/feassemble/GeometryLine3D.hh"

#include "data/QuadratureData1Din3DLinear.hh"
#include "data/QuadratureData1Din3DQuadratic.hh"

// ----------------------------------------------------------------------
CPPUNIT_TEST_SUITE_REGISTRATION( pylith::feassemble::TestQuadrature1Din3D );

// ----------------------------------------------------------------------
// Test constructor
void
pylith::feassemble::TestQuadrature1Din3D::testConstructor(void)
{ // testConstructor
  QuadratureRefCell refCell;
  Quadrature1Din3D quadrature(refCell);
} // testConstructor

// ----------------------------------------------------------------------
// Test computeGeometry() w/linear basis fns
void
pylith::feassemble::TestQuadrature1Din3D::testLinear(void)
{ // testLinear
  GeometryLine3D geometry;
  QuadratureRefCell refCell;
  refCell.refGeometry(&geometry);

  Quadrature1Din3D q(refCell);
  QuadratureData1Din3DLinear data;

  _testComputeGeometry(&q, &refCell, data);
} // testLinear

// ----------------------------------------------------------------------
// Test computeGeometry() w/quadratic basis fns
void
pylith::feassemble::TestQuadrature1Din3D::testQuadratic(void)
{ // testQuadratic
  GeometryLine3D geometry;
  QuadratureRefCell refCell;
  refCell.refGeometry(&geometry);

  Quadrature1Din3D q(refCell);
  QuadratureData1Din3DQuadratic data;

  _testComputeGeometry(&q, &refCell, data);
} // testQuadratic

// End of file 
