// -*- C++ -*-
//
// ----------------------------------------------------------------------
//
//                           Brad T. Aagaard
//                        U.S. Geological Survey
//
// {LicenseText}
//
// ----------------------------------------------------------------------
//

#include <portinfo>

#include "TestQuadrature.hh" // Implementation of class methods

#include "pylith/feassemble/Quadrature1D.hh" // USES Quadrature1D

#include <sstream> // USES std::stringstream

// ----------------------------------------------------------------------
CPPUNIT_TEST_SUITE_REGISTRATION( pylith::feassemble::TestQuadrature );

// ----------------------------------------------------------------------
// Test clone
void
pylith::feassemble::TestQuadrature::testClone(void)
{ // testClone
  Quadrature1D qOrig;
  
  const double jacobianTol = 1.0;

  const int cellDim = 1;
  const int numCorners = 2;
  const int numQuadPts = 1;
  const int spaceDim = 1;
  const double basis[] = { 0.5, 0.5 };
  const double basisDeriv[] = { -1.0, 1.0 };
  const double quadPtsRef[] = { 0.5 };
  const double quadWts[] = { 1.0 };

  // Semi-random values manually set to check cloning
  const double quadPts[] = { 2.0 };
  const double jacobian[] = { 4.0 };
  const double jacobianInv[] = { 0.25 };
  const double jacobianDet[] = { 5.0 };

  qOrig.jacobianTolerance(jacobianTol);
  qOrig.initialize(basis, basisDeriv, quadPtsRef, quadWts,
		   cellDim, numCorners, numQuadPts, spaceDim);
  int size = 1;
  memcpy(qOrig._quadPts, quadPts, size*sizeof(double));
  memcpy(qOrig._jacobian, jacobian, size*sizeof(double));
  memcpy(qOrig._jacobianInv, jacobianInv, size*sizeof(double));
  memcpy(qOrig._jacobianDet, jacobianDet, size*sizeof(double));
  
  const Quadrature* qCopy = qOrig.clone();
  CPPUNIT_ASSERT(0 != qCopy);

  // Make sure values match those passed
  CPPUNIT_ASSERT_EQUAL(jacobianTol, qCopy->_jacobianTol);

  CPPUNIT_ASSERT_EQUAL(cellDim, qCopy->_cellDim);
  CPPUNIT_ASSERT_EQUAL(numCorners, qCopy->_numCorners);
  CPPUNIT_ASSERT_EQUAL(numQuadPts, qCopy->_numQuadPts);
  CPPUNIT_ASSERT_EQUAL(spaceDim, qCopy->_spaceDim);

  CPPUNIT_ASSERT(0 != qCopy->_basis);
  size = numCorners * numQuadPts;
  for (int i=0; i < size; ++i)
    CPPUNIT_ASSERT_EQUAL(basis[i], qCopy->_basis[i]);

  CPPUNIT_ASSERT(0 != qCopy->_basisDeriv);
  size = numCorners * numQuadPts * spaceDim;
  for (int i=0; i < size; ++i)
    CPPUNIT_ASSERT_EQUAL(basisDeriv[i], qCopy->_basisDeriv[i]);

  CPPUNIT_ASSERT(0 != qCopy->_quadPtsRef);
  size = numQuadPts * cellDim;
  for (int i=0; i < size; ++i)
    CPPUNIT_ASSERT_EQUAL(quadPtsRef[i], qCopy->_quadPtsRef[i]);

  CPPUNIT_ASSERT(0 != qCopy->_quadWts);
  size = numQuadPts;
  for (int i=0; i < size; ++i)
    CPPUNIT_ASSERT_EQUAL(quadWts[i], qCopy->_quadWts[i]);

  size = 1;

  CPPUNIT_ASSERT(0 != qCopy->_quadPts);
  for (int i=0; i < size; ++i)
    CPPUNIT_ASSERT_EQUAL(quadPts[i], qCopy->_quadPts[i]);
  CPPUNIT_ASSERT(0 != qCopy->_jacobian);
  for (int i=0; i < size; ++i)
    CPPUNIT_ASSERT_EQUAL(jacobian[i], qCopy->_jacobian[i]);
  CPPUNIT_ASSERT(0 != qCopy->_jacobianInv);
  for (int i=0; i < size; ++i)
    CPPUNIT_ASSERT_EQUAL(jacobianInv[i], qCopy->_jacobianInv[i]);
  CPPUNIT_ASSERT(0 != qCopy->_jacobianDet);
  for (int i=0; i < size; ++i)
    CPPUNIT_ASSERT_EQUAL(jacobianDet[i], qCopy->_jacobianDet[i]);
} // testCopy

// ----------------------------------------------------------------------
// Test jacobianTolerance()
void
pylith::feassemble::TestQuadrature::testJacobianTol(void)
{ // testJacobianTol
  Quadrature1D q;
  const double tolerance = 1.0;
  q.jacobianTolerance(tolerance);
  CPPUNIT_ASSERT_EQUAL(tolerance, q._jacobianTol);
} // testJacobianTol

// ----------------------------------------------------------------------
// Test initialize()
void
pylith::feassemble::TestQuadrature::testInitialize(void)
{ // initialize
  Quadrature1D q;
  
  const int cellDim = 1;
  const int numCorners = 2;
  const int numQuadPts = 1;
  const int spaceDim = 1;
  const double basis[] = { 0.5, 0.5 };
  const double basisDeriv[] = { -1.0, 1.0 };
  const double quadPtsRef[] = { 0.5 };
  const double quadWts[] = { 1.0 };
  const double jacobianTol = 1.0;

  q.initialize(basis, basisDeriv, quadPtsRef, quadWts,
	       cellDim, numCorners, numQuadPts, spaceDim);
  
  CPPUNIT_ASSERT_EQUAL(cellDim, q._cellDim);
  CPPUNIT_ASSERT_EQUAL(numCorners, q._numCorners);
  CPPUNIT_ASSERT_EQUAL(numQuadPts, q._numQuadPts);
  CPPUNIT_ASSERT_EQUAL(spaceDim, q._spaceDim);

  int size = numCorners * numQuadPts;
  for (int i=0; i < size; ++i)
    CPPUNIT_ASSERT_EQUAL(basis[i], q._basis[i]);

  size = numCorners * numQuadPts * spaceDim;
  for (int i=0; i < size; ++i)
    CPPUNIT_ASSERT_EQUAL(basisDeriv[i], q._basisDeriv[i]);

  size = numQuadPts * cellDim;
  for (int i=0; i < size; ++i)
    CPPUNIT_ASSERT_EQUAL(quadPtsRef[i], q._quadPtsRef[i]);

  size = numQuadPts;
  for (int i=0; i < size; ++i)
    CPPUNIT_ASSERT_EQUAL(quadWts[i], q._quadWts[i]);

  // Make sure Jacobian stuff has been allocated
  CPPUNIT_ASSERT(0 != q._jacobian);
  CPPUNIT_ASSERT(0 != q._jacobianInv);
  CPPUNIT_ASSERT(0 != q._jacobianDet);
  CPPUNIT_ASSERT(0 != q._quadPts);
} // initialize

// End of file 
