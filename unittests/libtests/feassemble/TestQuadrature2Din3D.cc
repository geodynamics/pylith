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

#include "TestQuadrature2Din3D.hh" // Implementation of class methods

#include "pylith/feassemble/Quadrature2Din3D.hh"
#include "pylith/feassemble/QuadratureRefCell.hh"
#include "pylith/feassemble/GeometryTri3D.hh"

#include "data/QuadratureData2Din3DLinearXYZ.hh"
#include "data/QuadratureData2Din3DLinearXY.hh"
#include "data/QuadratureData2Din3DLinearYZ.hh"
#include "data/QuadratureData2Din3DLinearXZ.hh"
#include "data/QuadratureData2Din3DQuadratic.hh"

#include "pylith/utils/error.h" // USES PYLITH_METHOD_BEGIN/END

// ----------------------------------------------------------------------
CPPUNIT_TEST_SUITE_REGISTRATION( pylith::feassemble::TestQuadrature2Din3D );

// ----------------------------------------------------------------------
// Test constructor
void
pylith::feassemble::TestQuadrature2Din3D::testConstructor(void)
{ // testConstructor
  PYLITH_METHOD_BEGIN;

  QuadratureRefCell refCell;
  Quadrature2Din3D quadrature(refCell);

  PYLITH_METHOD_END;
} // testConstructor

// ----------------------------------------------------------------------
// Test computeGeometry() w/linear basis fns
void
pylith::feassemble::TestQuadrature2Din3D::testLinearXYZ(void)
{ // testLinearXYZ
  PYLITH_METHOD_BEGIN;

  GeometryTri3D geometry;
  QuadratureRefCell refCell;
  refCell.refGeometry(&geometry);

  Quadrature2Din3D q(refCell);
  QuadratureData2Din3DLinearXYZ data;

  _testComputeGeometry(&q, &refCell, data);

  PYLITH_METHOD_END;
} // testLinearXYZ

// ----------------------------------------------------------------------
// Test computeGeometry() w/linear basis fns
void
pylith::feassemble::TestQuadrature2Din3D::testLinearXY(void)
{ // testLinearXY
  PYLITH_METHOD_BEGIN;

  GeometryTri3D geometry;
  QuadratureRefCell refCell;
  refCell.refGeometry(&geometry);

  Quadrature2Din3D q(refCell);
  QuadratureData2Din3DLinearXY data;

  _testComputeGeometry(&q, &refCell, data);

  PYLITH_METHOD_END;
} // testLinearXY

// ----------------------------------------------------------------------
// Test computeGeometry() w/linear basis fns
void
pylith::feassemble::TestQuadrature2Din3D::testLinearYZ(void)
{ // testLinearYZ
  PYLITH_METHOD_BEGIN;

  GeometryTri3D geometry;
  QuadratureRefCell refCell;
  refCell.refGeometry(&geometry);

  Quadrature2Din3D q(refCell);
  QuadratureData2Din3DLinearYZ data;

  _testComputeGeometry(&q, &refCell, data);

  PYLITH_METHOD_END;
} // testLinearYZ

// ----------------------------------------------------------------------
// Test computeGeometry() w/linear basis fns
void
pylith::feassemble::TestQuadrature2Din3D::testLinearXZ(void)
{ // testLinearXZ
  PYLITH_METHOD_BEGIN;

  GeometryTri3D geometry;
  QuadratureRefCell refCell;
  refCell.refGeometry(&geometry);

  Quadrature2Din3D q(refCell);
  QuadratureData2Din3DLinearXZ data;

  _testComputeGeometry(&q, &refCell, data);

  PYLITH_METHOD_END;
} // testLinearXZ

// ----------------------------------------------------------------------
// Test computeGeometry() w/quadratic basis fns
void
pylith::feassemble::TestQuadrature2Din3D::testQuadratic(void)
{ // testQuadratic
  PYLITH_METHOD_BEGIN;

  GeometryTri3D geometry;
  QuadratureRefCell refCell;
  refCell.refGeometry(&geometry);

  Quadrature2Din3D q(refCell);
  QuadratureData2Din3DQuadratic data;

  _testComputeGeometry(&q, &refCell, data);

  PYLITH_METHOD_END;
} // testQuadratic


// End of file 
