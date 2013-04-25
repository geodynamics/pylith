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
// Copyright (c) 2010-2013 University of California, Davis
//
// See COPYING for license information.
//
// ----------------------------------------------------------------------
//

#include <portinfo>

#include "TestQuadrature0D.hh" // Implementation of class methods

#include "pylith/feassemble/Quadrature0D.hh"
#include "pylith/feassemble/QuadratureRefCell.hh"

#include "pylith/utils/error.h" // USES PYLITH_METHOD_BEGIN/END

// ----------------------------------------------------------------------
CPPUNIT_TEST_SUITE_REGISTRATION( pylith::feassemble::TestQuadrature0D );

// ----------------------------------------------------------------------
// Test constructor
void
pylith::feassemble::TestQuadrature0D::testConstructor(void)
{ // testConstructor
  PYLITH_METHOD_BEGIN;

  QuadratureRefCell refCell;
  Quadrature0D quadrature(refCell);

  PYLITH_METHOD_END;
} // testConstructor

// ----------------------------------------------------------------------
// Test computeGeometry().
void
pylith::feassemble::TestQuadrature0D::testPoint(void)
{ // testPoint
  PYLITH_METHOD_BEGIN;

  const int cellDim = 0;
  const int numBasis = 1;
  const int numQuadPts = 1;
  const int spaceDim = 1;
  const PylithScalar basis[] = { 1.0 };
  const PylithScalar basisDerivRef[] = { 1.0 };
  const PylithScalar quadPtsRef[] = { 0.0 };
  const PylithScalar quadWts[] = { 1.0 };

  const int numVertices = 1;
  const int numCells = 1;
  const PylithScalar vertCoordsData[] = { 1.1 };
  const int cells[] = { 0 };
  const PylithScalar quadPts[] = { 1.1 };
  const PylithScalar jacobian[] = { 1.0 };
  const PylithScalar jacobianInv[] = { 1.0 };
  const PylithScalar jacobianDet[] = { 1.0 };

  scalar_array vertCoords(vertCoordsData, numBasis*spaceDim);

  const PylithScalar minJacobian = 1.0e-06;
  
  QuadratureRefCell refCell;
  refCell.minJacobian(minJacobian);
  refCell.initialize(basis, numQuadPts, numBasis,
		     basisDerivRef, numQuadPts, numBasis, cellDim,
		     quadPtsRef, numQuadPts, cellDim,
		     quadWts, numQuadPts,
		     spaceDim);

  Quadrature0D engine(refCell);

  engine.initialize();
  engine.computeGeometry(vertCoords, 0);

  const PylithScalar tolerance = 1.0e-06;
  CPPUNIT_ASSERT_DOUBLES_EQUAL(quadPts[0], engine._quadPts[0], tolerance);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(jacobian[0], engine._jacobian[0], tolerance);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(jacobianInv[0], engine._jacobianInv[0], tolerance);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(jacobianDet[0], engine._jacobianDet[0], tolerance);

  PYLITH_METHOD_END;
} // testPoint


// End of file 
