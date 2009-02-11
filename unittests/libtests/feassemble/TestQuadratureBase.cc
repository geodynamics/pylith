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

#include "TestQuadratureBase.hh" // Implementation of class methods

#include "pylith/feassemble/QuadratureBase.hh" // USES QuadratureBase
#include "pylith/feassemble/GeometryLine1D.hh" // USES GeometryLine1D

#include "data/QuadratureData.hh" // USES QuadratureData

#include <string.h> // USES memcpy()

// ----------------------------------------------------------------------
CPPUNIT_TEST_SUITE_REGISTRATION( pylith::feassemble::TestQuadratureBase );

// ----------------------------------------------------------------------
// Test constructor.
void
pylith::feassemble::TestQuadratureBase::testConstructor(void)
{ // testConstructor
  QuadratureBase q;
} // testMinJacobian

// ----------------------------------------------------------------------
// Test minJacobian()
void
pylith::feassemble::TestQuadratureBase::testMinJacobian(void)
{ // testMinJacobian
  QuadratureBase q;

  const double min = 1.0;
  q.minJacobian(min);
  CPPUNIT_ASSERT_EQUAL(min, q._minJacobian);
} // testMinJacobian

// ----------------------------------------------------------------------
// Test refGeometry()
void
pylith::feassemble::TestQuadratureBase::testRefGeometry(void)
{ // testRefGeometry
  GeometryLine1D geometry;
  QuadratureBase quadrature;

  quadrature.refGeometry(&geometry);
  const CellGeometry& test = quadrature.refGeometry();

  CPPUNIT_ASSERT_EQUAL(geometry.cellDim(), test.cellDim());
  CPPUNIT_ASSERT_EQUAL(geometry.spaceDim(), test.spaceDim());
  CPPUNIT_ASSERT_EQUAL(geometry.numCorners(), test.numCorners());
} // testRefGeometry

// ----------------------------------------------------------------------
// Test initialize()
void
pylith::feassemble::TestQuadratureBase::testInitialize(void)
{ // initialize
  
  const int cellDim = 1;
  const int numBasis = 2;
  const int numQuadPts = 1;
  const int spaceDim = 1;
  const double basis[] = { 0.5, 0.5 };
  const double basisDerivRef[] = { -0.5, 0.5 };
  const double quadPtsRef[] = { 0.0 };
  const double quadWts[] = { 2.0 };
  const double minJacobian = 1.0;

  QuadratureBase q;
  q.initialize(basis, basisDerivRef, quadPtsRef, quadWts,
	       cellDim, numBasis, numQuadPts, spaceDim);
  
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
} // initialize


// End of file 
