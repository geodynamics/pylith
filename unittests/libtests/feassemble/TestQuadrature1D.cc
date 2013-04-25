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

#include "TestQuadrature1D.hh" // Implementation of class methods

#include "pylith/feassemble/Quadrature1D.hh"
#include "pylith/feassemble/QuadratureRefCell.hh"
#include "pylith/feassemble/GeometryLine1D.hh"

#include "data/QuadratureData1DLinear.hh"
#include "data/QuadratureData1DQuadratic.hh"

#include "pylith/utils/petscerror.h" // USES PYLITH_METHOD_BEGIN/END

// ----------------------------------------------------------------------
CPPUNIT_TEST_SUITE_REGISTRATION( pylith::feassemble::TestQuadrature1D );

// ----------------------------------------------------------------------
// Test constructor
void
pylith::feassemble::TestQuadrature1D::testConstructor(void)
{ // testConstructor
  PYLITH_METHOD_BEGIN;

  QuadratureRefCell refCell;
  Quadrature1D quadrature(refCell);

  PYLITH_METHOD_END;
} // testConstructor

// ----------------------------------------------------------------------
// Test computeGeometry() w/linear basis fns
void
pylith::feassemble::TestQuadrature1D::testLinear(void)
{ // testLinear
  PYLITH_METHOD_BEGIN;

  GeometryLine1D geometry;
  QuadratureRefCell refCell;
  refCell.refGeometry(&geometry);

  Quadrature1D q(refCell);
  QuadratureData1DLinear data;

  _testComputeGeometry(&q, &refCell, data);

  PYLITH_METHOD_END;
} // testLinear

// ----------------------------------------------------------------------
// Test computeGeometry() w/quadratic basis fns
void
pylith::feassemble::TestQuadrature1D::testQuadratic(void)
{ // testQuadratic
  PYLITH_METHOD_BEGIN;

  GeometryLine1D geometry;
  QuadratureRefCell refCell;
  refCell.refGeometry(&geometry);

  Quadrature1D q(refCell);
  QuadratureData1DQuadratic data;

  _testComputeGeometry(&q, &refCell, data);

  PYLITH_METHOD_END;
} // testQuadratic


// End of file 
