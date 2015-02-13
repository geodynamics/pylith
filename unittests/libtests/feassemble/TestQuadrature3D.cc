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
// Copyright (c) 2010-2015 University of California, Davis
//
// See COPYING for license information.
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

#include "pylith/utils/error.h" // USES PYLITH_METHOD_BEGIN/END

// ----------------------------------------------------------------------
CPPUNIT_TEST_SUITE_REGISTRATION( pylith::feassemble::TestQuadrature3D );

// ----------------------------------------------------------------------
// Test constructor
void
pylith::feassemble::TestQuadrature3D::testConstructor(void)
{ // testConstructor
  PYLITH_METHOD_BEGIN;

  QuadratureRefCell refCell;
  Quadrature3D quadrature(refCell);

  PYLITH_METHOD_END;
} // testConstructor

// ----------------------------------------------------------------------
// Test computeGeometry() w/linear basis fns
void
pylith::feassemble::TestQuadrature3D::testLinear(void)
{ // testLinear
  PYLITH_METHOD_BEGIN;

  GeometryTet3D geometry;
  QuadratureRefCell refCell;
  refCell.refGeometry(&geometry);

  Quadrature3D q(refCell);
  QuadratureData3DLinear data;

  _testComputeGeometry(&q, &refCell, data);

  PYLITH_METHOD_END;
} // testLinear

// ----------------------------------------------------------------------
// Test computeGeometry() w/quadratic basis fns
void
pylith::feassemble::TestQuadrature3D::testQuadratic(void)
{ // testQuadratic
  PYLITH_METHOD_BEGIN;

  GeometryTet3D geometry;
  QuadratureRefCell refCell;
  refCell.refGeometry(&geometry);

  Quadrature3D q(refCell);
  QuadratureData3DQuadratic data;

  _testComputeGeometry(&q, &refCell, data);

  PYLITH_METHOD_END;
} // testQuadratic


// End of file 
