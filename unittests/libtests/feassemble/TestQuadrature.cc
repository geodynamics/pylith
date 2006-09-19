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
  // Semi-random values manually set to check cloning
  const double jacobianTol = 1.0;
  const int cellDim = 1;
  const int numCorners = 2;
  const int numQuadPts = 1;
  const int spaceDim = 1;
  const double basis[] = { 0.2, 0.4 };
  const double basisDeriv[] = { 0.8, 1.6 };
  const double quadPtsRef[] = { 3.2 };
  const double quadWts[] = { 6.4 };
  const double quadPts[] = { 12.8 };
  const double jacobian[] = { 2.56 };
  const double jacobianInv[] = { 5.12 };
  const double jacobianDet[] = { 10.24 };

  // Set values
  Quadrature1D qOrig;
  qOrig._jacobianTol = jacobianTol;
  qOrig._cellDim = cellDim;
  qOrig._numCorners = numCorners;
  qOrig._numQuadPts = numQuadPts;
  qOrig._spaceDim = spaceDim;

  int size = 2;
  qOrig._basis = (size > 0) ? new double[size] : 0;
  memcpy(qOrig._basis, basis, size*sizeof(double));
  
  size = 2;
  qOrig._basisDeriv = (size > 0) ? new double[size] : 0;
  memcpy(qOrig._basisDeriv, basisDeriv, size*sizeof(double));

  size = 1;
  qOrig._quadPtsRef = (size > 0) ? new double[size] : 0;
  memcpy(qOrig._quadPtsRef, quadPtsRef, size*sizeof(double));

  size = 1;
  qOrig._quadWts = (size > 0) ? new double[size] : 0;
  memcpy(qOrig._quadWts, quadWts, size*sizeof(double));

  size = 1;
  qOrig._quadPts = (size > 0) ? new double[size] : 0;
  memcpy(qOrig._quadPts, quadPts, size*sizeof(double));

  size = 1;
  qOrig._jacobian = (size > 0) ? new double[size] : 0;
  memcpy(qOrig._jacobian, jacobian, size*sizeof(double));

  size = 1;
  qOrig._jacobianInv = (size > 0) ? new double[size] : 0;
  memcpy(qOrig._jacobianInv, jacobianInv, size*sizeof(double));

  size = 1;
  qOrig._jacobianDet = (size > 0) ? new double[size] : 0;
  memcpy(qOrig._jacobianDet, jacobianDet, size*sizeof(double));

  // Clone
  const Quadrature* qCopy = qOrig.clone();

  // Check clone
  CPPUNIT_ASSERT(0 != qCopy);

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
  
  const int cellDim = 1;
  const int numCorners = 2;
  const int numQuadPts = 1;
  const int spaceDim = 1;
  const double basis[] = { 0.5, 0.5 };
  const double basisDeriv[] = { -1.0, 1.0 };
  const double quadPtsRef[] = { 0.5 };
  const double quadWts[] = { 1.0 };
  const double jacobianTol = 1.0;

  Quadrature1D q;
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
