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

#include "TestQuadratureEngine.hh" // Implementation of class methods

#include "pylith/feassemble/QuadratureRefCell.hh" // USES QuadratureRefCell
#include "pylith/feassemble/Quadrature1D.hh" // USES Quadrature1D

#include "data/QuadratureData.hh" // USES QuadratureData

#include <string.h> // USES memcpy()

// ----------------------------------------------------------------------
CPPUNIT_TEST_SUITE_REGISTRATION( pylith::feassemble::TestQuadratureEngine );

// ----------------------------------------------------------------------
// Test copy constructor.
void
pylith::feassemble::TestQuadratureEngine::testCopyConstructor(void)
{ // testClone
  // Semi-random values manually set to check cloning
  const double quadPtsE[] = { 12.8 };
  const double jacobianE[] = { 2.56 };
  const double jacobianInvE[] = { 5.12 };
  const double jacobianDetE[] = { 10.24 };
  const double basisDerivE[] = { 0.8, 1.6 };

  QuadratureRefCell refCell;

  Quadrature1D engineOrig(refCell);
  size_t size = 0;

  // Set values
  size = 1;
  engineOrig._quadPts.resize(size);
  memcpy(&engineOrig._quadPts[0], quadPtsE, size*sizeof(double));

  size = 1;
  engineOrig._jacobian.resize(size);
  memcpy(&engineOrig._jacobian[0], jacobianE, size*sizeof(double));

  size = 1;
  engineOrig._jacobianInv.resize(size);
  memcpy(&engineOrig._jacobianInv[0], jacobianInvE, size*sizeof(double));

  size = 1;
  engineOrig._jacobianDet.resize(size);
  memcpy(&engineOrig._jacobianDet[0], jacobianDetE, size*sizeof(double));

  // Copy
  const QuadratureEngine* engine = engineOrig.clone();
  CPPUNIT_ASSERT(0 != engine);

  const double_array& quadPts = engine->quadPts();
  size = 1;
  CPPUNIT_ASSERT_EQUAL(size, quadPts.size());
  for (int i=0; i < size; ++i)
    CPPUNIT_ASSERT_EQUAL(quadPtsE[i], quadPts[i]);

  const double_array& jacobian = engine->_jacobian;
  size = 1;
  CPPUNIT_ASSERT_EQUAL(size, jacobian.size());
  for (int i=0; i < size; ++i)
    CPPUNIT_ASSERT_EQUAL(jacobianE[i], jacobian[i]);

  const double_array& jacobianInv = engine->_jacobianInv;
  size = 1;
  CPPUNIT_ASSERT_EQUAL(size, jacobianInv.size());
  for (int i=0; i < size; ++i)
    CPPUNIT_ASSERT_EQUAL(jacobianInvE[i], jacobianInv[i]);

  const double_array& jacobianDet = engine->jacobianDet();
  size = 1;
  CPPUNIT_ASSERT_EQUAL(size, jacobianDet.size());
  for (int i=0; i < size; ++i)
    CPPUNIT_ASSERT_EQUAL(jacobianDetE[i], jacobianDet[i]);

  delete engine; engine = 0;
} // testCopyConstructor

// ----------------------------------------------------------------------
// Test initialize().
void
pylith::feassemble::TestQuadratureEngine::testInitialize(void)
{ // testInitialize
  const int cellDim = 2;
  const int numBasis = 5;
  const int numQuadPts = 1;
  const int spaceDim = 3;
  const double basis[] = { 
    1.1, 1.2, 1.3, 1.4, 1.5
  };
  const double basisDerivRef[] = {
    2.1, 2.2, 2.3,
    2.4, 2.5, 2.6,
    2.7, 2.8, 2.9,
    2.10, 2.11, 2.12,
    2.13, 2.14, 2.15,
  };
  const double quadPtsRef[] = { 3.1, 3.2, 3.3 };
  const double quadWts[] = { 4.0 };

  QuadratureRefCell refCell;
  refCell.initialize(basis, numQuadPts, numBasis,
		     basisDerivRef, numQuadPts, numBasis, cellDim,
		     quadPtsRef, numQuadPts, cellDim,
		     quadWts, numQuadPts,
		     spaceDim);

  Quadrature1D engine(refCell);
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
} // testInitialize

// ----------------------------------------------------------------------
// Test computeGeometry().
void
pylith::feassemble::TestQuadratureEngine::_testComputeGeometry(
					      QuadratureEngine* engine,
					      QuadratureRefCell* refCell,
					      const QuadratureData& data) const
{ // testComputeGeometry
  const int cellDim = data.cellDim;
  const int numBasis = data.numBasis;
  const int numQuadPts = data.numQuadPts;
  const int spaceDim = data.spaceDim;
  const double* basis = data.basis;
  const double* basisDerivRef = data.basisDerivRef;
  const double* quadPtsRef = data.quadPtsRef;
  const double* quadWts = data.quadWts;

  const int numVertices = data.numVertices;
  const int numCells = data.numCells;
  const double_array vertCoords(data.vertices, numBasis*spaceDim);
  const int* cells = data.cells;
  const double* quadPtsE = data.quadPts;
  const double* jacobianE = data.jacobian;
  const double* jacobianInvE = data.jacobianInv;
  const double* jacobianDetE = data.jacobianDet;
  const double* basisDerivE = data.basisDeriv;

  CPPUNIT_ASSERT(1 == numCells);
  CPPUNIT_ASSERT(0 != engine);
  CPPUNIT_ASSERT(0 != refCell);

  const double minJacobian = 1.0e-06;
  refCell->minJacobian(minJacobian);
  refCell->initialize(basis, numQuadPts, numBasis,
		      basisDerivRef, numQuadPts, numBasis, cellDim,
		      quadPtsRef, numQuadPts, cellDim,
		      quadWts, numQuadPts,
		      spaceDim);

  // Check values from computeGeometry()
  engine->initialize();
  engine->computeGeometry(vertCoords, 0);

  const double tolerance = 1.0e-06;
  int size = numQuadPts * spaceDim;
  for (int i=0; i < size; ++i)
    CPPUNIT_ASSERT_DOUBLES_EQUAL(quadPtsE[i], engine->_quadPts[i], tolerance);

  size = numQuadPts * cellDim * spaceDim;
  for (int i=0; i < size; ++i)
    CPPUNIT_ASSERT_DOUBLES_EQUAL(jacobianE[i], engine->_jacobian[i], 
				 tolerance);

  size = numQuadPts * spaceDim * cellDim;
  for (int i=0; i < size; ++i)
    CPPUNIT_ASSERT_DOUBLES_EQUAL(jacobianInvE[i], engine->_jacobianInv[i], 
				 tolerance);

  size = numQuadPts;
  for (int i=0; i < size; ++i)
    CPPUNIT_ASSERT_DOUBLES_EQUAL(jacobianDetE[i], engine->_jacobianDet[i], 
				 tolerance);

  size = numQuadPts * numBasis * spaceDim;
  for (int i=0; i < size; ++i)
    CPPUNIT_ASSERT_DOUBLES_EQUAL(basisDerivE[i], engine->_basisDeriv[i], 
				 tolerance);
} // testComputeGeometry


// End of file 
