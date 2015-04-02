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

#include "TestQuadratureEngine.hh" // Implementation of class methods

#include "pylith/feassemble/QuadratureRefCell.hh" // USES QuadratureRefCell
#include "pylith/feassemble/Quadrature2D.hh" // USES Quadrature2D

#include "data/QuadratureData.hh" // USES QuadratureData

#include "pylith/utils/error.h" // USES PYLITH_METHOD_BEGIN/END

#include <string.h> // USES memcpy()

// ----------------------------------------------------------------------
CPPUNIT_TEST_SUITE_REGISTRATION( pylith::feassemble::TestQuadratureEngine );

// ----------------------------------------------------------------------
// Test copy constructor.
void
pylith::feassemble::TestQuadratureEngine::testCopyConstructor(void)
{ // testClone
  PYLITH_METHOD_BEGIN;

  // Semi-random values manually set to check cloning
  const PylithScalar quadPtsE[] = { 12.8 };
  const PylithScalar jacobianE[] = { 2.56 };
  const PylithScalar jacobianInvE[] = { 5.12 };
  const PylithScalar jacobianDetE[] = { 10.24 };
  const PylithScalar basisDerivE[] = { 0.8, 1.6 };

  QuadratureRefCell refCell;

  Quadrature2D engineOrig(refCell);
  size_t size = 0;

  // Set values
  size = 1;
  engineOrig._quadPts.resize(size);
  memcpy(&engineOrig._quadPts[0], quadPtsE, size*sizeof(PylithScalar));

  size = 1;
  engineOrig._jacobian.resize(size);
  memcpy(&engineOrig._jacobian[0], jacobianE, size*sizeof(PylithScalar));

  size = 1;
  engineOrig._jacobianInv.resize(size);
  memcpy(&engineOrig._jacobianInv[0], jacobianInvE, size*sizeof(PylithScalar));

  size = 1;
  engineOrig._jacobianDet.resize(size);
  memcpy(&engineOrig._jacobianDet[0], jacobianDetE, size*sizeof(PylithScalar));

  // Copy
  const QuadratureEngine* engine = engineOrig.clone();
  CPPUNIT_ASSERT(engine);

  const scalar_array& quadPts = engine->quadPts();
  size = 1;
  CPPUNIT_ASSERT_EQUAL(size, quadPts.size());
  for (int i=0; i < size; ++i)
    CPPUNIT_ASSERT_EQUAL(quadPtsE[i], quadPts[i]);

  const scalar_array& jacobian = engine->_jacobian;
  size = 1;
  CPPUNIT_ASSERT_EQUAL(size, jacobian.size());
  for (int i=0; i < size; ++i)
    CPPUNIT_ASSERT_EQUAL(jacobianE[i], jacobian[i]);

  const scalar_array& jacobianInv = engine->_jacobianInv;
  size = 1;
  CPPUNIT_ASSERT_EQUAL(size, jacobianInv.size());
  for (int i=0; i < size; ++i)
    CPPUNIT_ASSERT_EQUAL(jacobianInvE[i], jacobianInv[i]);

  const scalar_array& jacobianDet = engine->jacobianDet();
  size = 1;
  CPPUNIT_ASSERT_EQUAL(size, jacobianDet.size());
  for (int i=0; i < size; ++i)
    CPPUNIT_ASSERT_EQUAL(jacobianDetE[i], jacobianDet[i]);

  delete engine; engine = 0;

  PYLITH_METHOD_END;
} // testCopyConstructor

// ----------------------------------------------------------------------
// Test initialize().
void
pylith::feassemble::TestQuadratureEngine::testInitialize(void)
{ // testInitialize
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

  size_t size = 0;

  size = numQuadPts * spaceDim;
  CPPUNIT_ASSERT_EQUAL(size, engine.quadPts().size());

  size = numQuadPts * cellDim * spaceDim;
  CPPUNIT_ASSERT_EQUAL(size, engine.jacobian().size());

  size = numQuadPts * cellDim * spaceDim;
  CPPUNIT_ASSERT_EQUAL(size, engine._jacobianInv.size());

  size = numQuadPts;
  CPPUNIT_ASSERT_EQUAL(size, engine.jacobianDet().size());

  size = numQuadPts * numBasis * spaceDim;
  CPPUNIT_ASSERT_EQUAL(size, engine.basisDeriv().size());

  PYLITH_METHOD_END;
} // testInitialize

// ----------------------------------------------------------------------
// Test computeGeometry().
void
pylith::feassemble::TestQuadratureEngine::_testComputeGeometry(
					      QuadratureEngine* engine,
					      QuadratureRefCell* refCell,
					      const QuadratureData& data) const
{ // testComputeGeometry
  PYLITH_METHOD_BEGIN;

  const int cellDim = data.cellDim;
  const int numBasis = data.numBasis;
  const int numQuadPts = data.numQuadPts;
  const int spaceDim = data.spaceDim;
  const PylithScalar* basis = data.basis;
  const PylithScalar* basisDerivRef = data.basisDerivRef;
  const PylithScalar* quadPtsRef = data.quadPtsRef;
  const PylithScalar* quadWts = data.quadWts;

  const int numVertices = data.numVertices;
  const int numCells = data.numCells;
  const scalar_array vertCoords(data.vertices, numBasis*spaceDim);
  const int* cells = data.cells;
  const PylithScalar* quadPtsE = data.quadPts;
  const PylithScalar* jacobianE = data.jacobian;
  const PylithScalar* jacobianInvE = data.jacobianInv;
  const PylithScalar* jacobianDetE = data.jacobianDet;
  const PylithScalar* basisDerivE = data.basisDeriv;

  CPPUNIT_ASSERT(1 == numCells);
  CPPUNIT_ASSERT(engine);
  CPPUNIT_ASSERT(refCell);

  const PylithScalar minJacobian = 1.0e-06;
  refCell->minJacobian(minJacobian);
  refCell->initialize(basis, numQuadPts, numBasis,
		      basisDerivRef, numQuadPts, numBasis, cellDim,
		      quadPtsRef, numQuadPts, cellDim,
		      quadWts, numQuadPts,
		      spaceDim);

  // Check values from computeGeometry()
  engine->initialize();
  engine->computeGeometry(&vertCoords[0], vertCoords.size(), 0);

  const PylithScalar tolerance = 1.0e-06;
  int size = numQuadPts * spaceDim;
  for (int i=0; i < size; ++i)
    CPPUNIT_ASSERT_DOUBLES_EQUAL(quadPtsE[i], engine->_quadPts[i], tolerance);

  size = numQuadPts * cellDim * spaceDim;
  for (int i=0; i < size; ++i)
    CPPUNIT_ASSERT_DOUBLES_EQUAL(jacobianE[i], engine->_jacobian[i], tolerance);

  size = numQuadPts * spaceDim * cellDim;
  for (int i=0; i < size; ++i)
    CPPUNIT_ASSERT_DOUBLES_EQUAL(jacobianInvE[i], engine->_jacobianInv[i], tolerance);

  size = numQuadPts;
  for (int i=0; i < size; ++i)
    CPPUNIT_ASSERT_DOUBLES_EQUAL(jacobianDetE[i], engine->_jacobianDet[i], tolerance);

  size = numQuadPts * numBasis * spaceDim;
  for (int i=0; i < size; ++i)
    CPPUNIT_ASSERT_DOUBLES_EQUAL(basisDerivE[i], engine->_basisDeriv[i], tolerance);

  PYLITH_METHOD_END;
} // testComputeGeometry


// End of file 
