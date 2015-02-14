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

#include "TestQuadratureRefCell.hh" // Implementation of class methods

#include "pylith/feassemble/QuadratureRefCell.hh" // USES QuadratureRefCell
#include "pylith/feassemble/GeometryLine2D.hh" // USES GeometryLine2D

#include "data/QuadratureData.hh" // USES QuadratureData

#include "pylith/utils/error.h" // USES PYLITH_METHOD_BEGIN/END

#include <string.h> // USES memcpy()

// ----------------------------------------------------------------------
CPPUNIT_TEST_SUITE_REGISTRATION( pylith::feassemble::TestQuadratureRefCell );

// ----------------------------------------------------------------------
// Test constructor.
void
pylith::feassemble::TestQuadratureRefCell::testConstructor(void)
{ // testConstructor
  PYLITH_METHOD_BEGIN;

  QuadratureRefCell q;

  PYLITH_METHOD_END;
} // testMinJacobian

// ----------------------------------------------------------------------
// Test minJacobian()
void
pylith::feassemble::TestQuadratureRefCell::testMinJacobian(void)
{ // testMinJacobian
  PYLITH_METHOD_BEGIN;

  QuadratureRefCell q;

  const PylithScalar min = 1.0;
  q.minJacobian(min);
  CPPUNIT_ASSERT_EQUAL(min, q._minJacobian);

  PYLITH_METHOD_END;
} // testMinJacobian

// ----------------------------------------------------------------------
// Test refGeometry()
void
pylith::feassemble::TestQuadratureRefCell::testRefGeometry(void)
{ // testRefGeometry
  PYLITH_METHOD_BEGIN;

  GeometryLine2D geometry;
  QuadratureRefCell quadrature;

  quadrature.refGeometry(&geometry);
  const CellGeometry& test = quadrature.refGeometry();

  CPPUNIT_ASSERT_EQUAL(geometry.cellDim(), test.cellDim());
  CPPUNIT_ASSERT_EQUAL(geometry.spaceDim(), test.spaceDim());
  CPPUNIT_ASSERT_EQUAL(geometry.numCorners(), test.numCorners());

  PYLITH_METHOD_END;
} // testRefGeometry

// ----------------------------------------------------------------------
// Test initialize()
void
pylith::feassemble::TestQuadratureRefCell::testInitialize(void)
{ // initialize
  PYLITH_METHOD_BEGIN;
  
  const int cellDim = 1;
  const int numBasis = 2;
  const int numQuadPts = 1;
  const int spaceDim = 1;
  const PylithScalar basis[] = { 0.5, 0.5 };
  const PylithScalar basisDerivRef[] = { -0.5, 0.5 };
  const PylithScalar quadPtsRef[] = { 0.0 };
  const PylithScalar quadWts[] = { 2.0 };
  const PylithScalar minJacobian = 1.0;

  QuadratureRefCell q;
  q.initialize(basis, numQuadPts, numBasis,
	       basisDerivRef, numQuadPts, numBasis, cellDim,
	       quadPtsRef, numQuadPts, cellDim,
	       quadWts, numQuadPts,
	       spaceDim);
  
  CPPUNIT_ASSERT_EQUAL(cellDim, q._cellDim);
  CPPUNIT_ASSERT_EQUAL(numBasis, q._numBasis);
  CPPUNIT_ASSERT_EQUAL(numQuadPts, q._numQuadPts);
  CPPUNIT_ASSERT_EQUAL(spaceDim, q._spaceDim);

  size_t size = numBasis * numQuadPts;
  for (int i=0; i < size; ++i)
    CPPUNIT_ASSERT_EQUAL(basis[i], q._basis[i]);

  size = numBasis * numQuadPts * spaceDim;
  for (int i=0; i < size; ++i)
    CPPUNIT_ASSERT_EQUAL(basisDerivRef[i], q._basisDerivRef[i]);

  size = numQuadPts * cellDim;
  for (int i=0; i < size; ++i)
    CPPUNIT_ASSERT_EQUAL(quadPtsRef[i], q._quadPtsRef[i]);

  size = numQuadPts;
  for (int i=0; i < size; ++i)
    CPPUNIT_ASSERT_EQUAL(quadWts[i], q._quadWts[i]);

  PYLITH_METHOD_END;
} // initialize


// End of file 
