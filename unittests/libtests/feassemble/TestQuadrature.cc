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

#include "TestQuadrature.hh" // Implementation of class methods

#include "pylith/topology/Mesh.hh" // USES Mesh
#include "pylith/feassemble/Quadrature.hh" // USES Quadrature
#include "pylith/feassemble/Quadrature2D.hh" // USES Quadrature1D

#include "pylith/feassemble/GeometryTri2D.hh" // USES GeometryTri2D

#include "data/QuadratureData2DLinear.hh" // USES QuadratureData2DLinear

#include <string.h> // USES memcpy()

// ----------------------------------------------------------------------
CPPUNIT_TEST_SUITE_REGISTRATION( pylith::feassemble::TestQuadrature );

// ----------------------------------------------------------------------
// Test copy constuctor.
void
pylith::feassemble::TestQuadrature::testCopyConstructor(void)
{ // testClone
  PYLITH_METHOD_BEGIN;

  // Semi-random values manually set to check cloning
  const PylithScalar minJacobianE = 1.0;
  const bool checkConditioning = true;
  const int cellDimE = 2;
  const int numBasisE = 3;
  const int numQuadPtsE = 1;
  const int spaceDimE = 2;
  const PylithScalar basisE[numQuadPtsE*numBasisE] = { 0.2, 0.3, 0.4 };
  const PylithScalar basisDerivE[numQuadPtsE*numBasisE*cellDimE] = { 0.8, 0.9, 1.2, 1.3, 1.6, 1.7 };
  const PylithScalar quadPtsRefE[numQuadPtsE*cellDimE] = { 3.2, 3.3 };
  const PylithScalar quadWtsE[numQuadPtsE] = { 6.4 };
  const PylithScalar quadPtsE[] = { 12.8 };
  const PylithScalar jacobianE[] = { 2.56 };
  const PylithScalar jacobianInvE[] = { 5.12 };
  const PylithScalar jacobianDetE[] = { 10.24 };
  GeometryTri2D geometry;

  // Set values
  Quadrature qOrig;
  qOrig.refGeometry(&geometry);
  qOrig.minJacobian(minJacobianE);
  qOrig.checkConditioning(checkConditioning);
  qOrig.initialize(basisE, numQuadPtsE, numBasisE,
		   basisDerivE, numQuadPtsE, numBasisE, cellDimE,
		   quadPtsRefE, numQuadPtsE, cellDimE,
		   quadWtsE, numQuadPtsE,
		   spaceDimE);

  // Copy
  Quadrature qCopy(qOrig);
  
  // Check copy
  CPPUNIT_ASSERT(!qCopy._engine);
  CPPUNIT_ASSERT_EQUAL(minJacobianE, qCopy._minJacobian);
  CPPUNIT_ASSERT_EQUAL(checkConditioning, qCopy._checkConditioning);
  CPPUNIT_ASSERT_EQUAL(cellDimE, qCopy.cellDim());
  CPPUNIT_ASSERT_EQUAL(numBasisE, qCopy.numBasis());
  CPPUNIT_ASSERT_EQUAL(numQuadPtsE, qCopy.numQuadPts());
  CPPUNIT_ASSERT_EQUAL(spaceDimE, qCopy.spaceDim());

  const scalar_array& basis = qCopy.basis();
  size_t size = numBasisE * numQuadPtsE;
  CPPUNIT_ASSERT_EQUAL(size, basis.size());
  for (int i=0; i < size; ++i)
    CPPUNIT_ASSERT_EQUAL(basisE[i], basis[i]);

  const scalar_array& basisDerivRef = qCopy._basisDerivRef;
  size = numBasisE * numQuadPtsE * spaceDimE;
  CPPUNIT_ASSERT_EQUAL(size, basisDerivRef.size());
  for (int i=0; i < size; ++i)
    CPPUNIT_ASSERT_EQUAL(basisDerivE[i], basisDerivRef[i]);

  const scalar_array& quadPtsRef = qCopy._quadPtsRef;
  size = numQuadPtsE * cellDimE;
  CPPUNIT_ASSERT_EQUAL(size, quadPtsRef.size());
  for (int i=0; i < size; ++i)
    CPPUNIT_ASSERT_EQUAL(quadPtsRefE[i], quadPtsRef[i]);

  const scalar_array& quadWts = qCopy.quadWts();
  size = numQuadPtsE;
  CPPUNIT_ASSERT_EQUAL(size, quadWts.size());
  for (int i=0; i < size; ++i)
    CPPUNIT_ASSERT_EQUAL(quadWtsE[i], quadWts[i]);

  CPPUNIT_ASSERT_EQUAL(geometry.cellDim(), qCopy.refGeometry().cellDim());
  CPPUNIT_ASSERT_EQUAL(geometry.spaceDim(), qCopy.refGeometry().spaceDim());
  CPPUNIT_ASSERT_EQUAL(geometry.numCorners(), qCopy.refGeometry().numCorners());

  PYLITH_METHOD_END;
} // testCopyConstructor

// ----------------------------------------------------------------------
// Test checkConditioning()
void
pylith::feassemble::TestQuadrature::testCheckConditioning(void)
{ // testCheckConditioning
  PYLITH_METHOD_BEGIN;

  Quadrature q;

  CPPUNIT_ASSERT_EQUAL(false, q.checkConditioning());
  q.checkConditioning(true);
  CPPUNIT_ASSERT_EQUAL(true, q.checkConditioning());
  q.checkConditioning(false);
  CPPUNIT_ASSERT_EQUAL(false, q.checkConditioning());

  PYLITH_METHOD_END;
} // testCheckConditioning

// ----------------------------------------------------------------------
// Test quadPts(), basisDeriv(), jacobian(), and jacobianDet().
void
pylith::feassemble::TestQuadrature::testEngineAccessors(void)
{ // testEngineAccessors
  PYLITH_METHOD_BEGIN;

  const int cellDim = 2;
  const int numBasis = 5;
  const int numQuadPts = 1;
  const int spaceDim = 3;
  const PylithScalar basis[] = { 
    1.1, 1.2, 1.3, 1.4, 1.5
  };
  const PylithScalar basisDerivRef[] = {
    2.1, 2.2, 2.3,
    2.4, 2.5, 2.6,
    2.7, 2.8, 2.9,
    2.10, 2.11, 2.12,
    2.13, 2.14, 2.15,
  };
  const PylithScalar quadPtsRef[] = { 3.1, 3.2, 3.3 };
  const PylithScalar quadWts[] = { 4.0 };

  QuadratureRefCell refCell;
  refCell.initialize(basis, numQuadPts, numBasis,
		     basisDerivRef, numQuadPts, numBasis, cellDim,
		     quadPtsRef, numQuadPts, cellDim,
		     quadWts, numQuadPts,
		     spaceDim);

  Quadrature2D engine(refCell);
  engine.initialize();

  Quadrature q;
  q._engine = engine.clone();

  size_t size = 0;

  size = numQuadPts * spaceDim;
  CPPUNIT_ASSERT_EQUAL(size, q.quadPts().size());

  size = numQuadPts * cellDim * spaceDim;
  CPPUNIT_ASSERT_EQUAL(size, q.jacobian().size());

  size = numQuadPts;
  CPPUNIT_ASSERT_EQUAL(size, q.jacobianDet().size());

  size = numQuadPts * numBasis * spaceDim;
  CPPUNIT_ASSERT_EQUAL(size, q.basisDeriv().size());

  PYLITH_METHOD_END;
} // testEngineAccessors

// ----------------------------------------------------------------------
// Test computeGeometry() will coordinates and cell.
void
pylith::feassemble::TestQuadrature::testComputeGeometryCell(void)
{ // testComputeGeometryCell
  PYLITH_METHOD_BEGIN;

  QuadratureData2DLinear data;
  const int cellDim = data.cellDim;
  const int numBasis = data.numBasis;
  const int numQuadPts = data.numQuadPts;
  const int spaceDim = data.spaceDim;

  const int numCells = data.numCells;
  const PylithScalar* vertCoords = data.vertices;
  const int vertCoordsSize = numBasis*spaceDim;
  const PylithScalar* quadPtsE = data.quadPts;
  const PylithScalar* jacobianE = data.jacobian;
  const PylithScalar* jacobianDetE = data.jacobianDet;
  const PylithScalar* basisDerivE = data.basisDeriv;

  const PylithScalar minJacobian = 1.0e-06;

  // Setup quadrature and compute geometry
  GeometryTri2D geometry;
  Quadrature quadrature;
  quadrature.refGeometry(&geometry);
  quadrature.minJacobian(minJacobian);
  quadrature.initialize(data.basis, numQuadPts, numBasis,
			data.basisDerivRef, numQuadPts, numBasis, cellDim,
			data.quadPtsRef, numQuadPts, cellDim,
			data.quadWts, numQuadPts,
			spaceDim);

  quadrature.initializeGeometry();
  quadrature.computeGeometry(vertCoords, vertCoordsSize, 0);
  
  size_t size = 0;

  // Check values from computeGeometry()
  const PylithScalar tolerance = 1.0e-06;

  const scalar_array& quadPts = quadrature.quadPts();
  size = numQuadPts * spaceDim;
  CPPUNIT_ASSERT_EQUAL(size, quadPts.size());
  for (size_t i=0; i < size; ++i)
    CPPUNIT_ASSERT_DOUBLES_EQUAL(quadPtsE[i], quadPts[i], tolerance);
  
  const scalar_array& jacobian = quadrature.jacobian();
  size = numQuadPts * cellDim * spaceDim;
  CPPUNIT_ASSERT_EQUAL(size, jacobian.size());
  for (size_t i=0; i < size; ++i)
    CPPUNIT_ASSERT_DOUBLES_EQUAL(jacobianE[i], jacobian[i], tolerance);
  
  const scalar_array& jacobianDet = quadrature.jacobianDet();
  size = numQuadPts;
  CPPUNIT_ASSERT_EQUAL(size, jacobianDet.size());
  for (size_t i=0; i < size; ++i)
    CPPUNIT_ASSERT_DOUBLES_EQUAL(jacobianDetE[i], jacobianDet[i], tolerance);
  
  const scalar_array& basisDeriv = quadrature.basisDeriv();
  size = numQuadPts * numBasis * spaceDim;
  CPPUNIT_ASSERT_EQUAL(size, basisDeriv.size());
  for (size_t i=0; i < size; ++i)
    CPPUNIT_ASSERT_DOUBLES_EQUAL(basisDerivE[i], basisDeriv[i], tolerance);

  quadrature.clear();

  CPPUNIT_ASSERT(!quadrature._engine);

  PYLITH_METHOD_END;
} // testComputeGeometryCell


// End of file 
