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

#include "TestIntegrator.hh" // Implementation of class methods

#include "pylith/feassemble/ElasticityExplicit.hh" // USES ElasticityExplicit
#include "pylith/feassemble/Quadrature1D.hh" // USES Quadrature1D

#include <petscmat.h>

#include <stdexcept> // TEMPORARY
// ----------------------------------------------------------------------
CPPUNIT_TEST_SUITE_REGISTRATION( pylith::feassemble::TestIntegrator );

// ----------------------------------------------------------------------
// Test quadrature().
void
pylith::feassemble::TestIntegrator::testQuadrature(void)
{ // testQuadrature
  // Since quadrature is cloned, test setting quadrature by testing
  // value of minJacobian

  Quadrature1D quadrature;
  const double minJacobian = 4.0;
  quadrature.minJacobian(minJacobian);
  
  ElasticityExplicit integrator;
  integrator.quadrature(&quadrature);

  CPPUNIT_ASSERT_EQUAL(minJacobian, integrator._quadrature->minJacobian());
} // testQuadrature

// ----------------------------------------------------------------------
// Test needNewJacobian().
void
pylith::feassemble::TestIntegrator::testNeedNewJacobian(void)
{ // testNeedNewJacobian
  ElasticityExplicit integrator;
  
  // Default should be false
  CPPUNIT_ASSERT_EQUAL(false, integrator._needNewJacobian);

  integrator._needNewJacobian = true;
  CPPUNIT_ASSERT_EQUAL(true, integrator._needNewJacobian);

  integrator._needNewJacobian = false;
  CPPUNIT_ASSERT_EQUAL(false, integrator._needNewJacobian);
} // testNeedNewJacobian

// ----------------------------------------------------------------------
// Test _initCellVector()
void
pylith::feassemble::TestIntegrator::testInitCellVector(void)
{ // testInitCellVector
  Quadrature1D quadrature;
  _initQuadrature(&quadrature);

  ElasticityExplicit integrator;
  integrator.quadrature(&quadrature);

  integrator._initCellVector();
  
  CPPUNIT_ASSERT(0 != integrator._cellVector);
  const int size = 
    quadrature.spaceDim() * quadrature.numBasis();
  for (int i=0; i < size; ++i)
    CPPUNIT_ASSERT_EQUAL(0.0, integrator._cellVector[i]);
} // testInitCellVector

// ----------------------------------------------------------------------
// Test _resetCellVector()
void
pylith::feassemble::TestIntegrator::testResetCellVector(void)
{ // testResetCellVector
  Quadrature1D quadrature;
  _initQuadrature(&quadrature);

  ElasticityExplicit integrator;
  integrator.quadrature(&quadrature);

  integrator._initCellVector();
  
  CPPUNIT_ASSERT(0 != integrator._cellVector);
  const int size = 
    quadrature.spaceDim() * quadrature.numBasis();
  for (int i=0; i < size; ++i)
    integrator._cellVector[i] = 1.4+2*i;
  integrator._resetCellVector();
  for (int i=0; i < size; ++i)
    CPPUNIT_ASSERT_EQUAL(0.0, integrator._cellVector[i]);
} // testResetCellVector

// ----------------------------------------------------------------------
// Test _initCellMatrix()
void
pylith::feassemble::TestIntegrator::testInitCellMatrix(void)
{ // testInitCellMatrix
  Quadrature1D quadrature;
  _initQuadrature(&quadrature);

  ElasticityExplicit integrator;
  integrator.quadrature(&quadrature);

  integrator._initCellMatrix();
  
  CPPUNIT_ASSERT(0 != integrator._cellMatrix);
  const int size = 
    quadrature.spaceDim() * quadrature.numBasis() *
    quadrature.spaceDim() * quadrature.numBasis();
  for (int i=0; i < size; ++i)
    CPPUNIT_ASSERT_EQUAL(0.0, integrator._cellMatrix[i]);
} // testInitCellMatrix

// ----------------------------------------------------------------------
// Test _resetCellMatrix()
void
pylith::feassemble::TestIntegrator::testResetCellMatrix(void)
{ // testResetCellMatrix
  Quadrature1D quadrature;
  _initQuadrature(&quadrature);

  ElasticityExplicit integrator;
  integrator.quadrature(&quadrature);

  integrator._initCellMatrix();
  
  CPPUNIT_ASSERT(0 != integrator._cellMatrix);
  const int size = 
    quadrature.spaceDim() * quadrature.numBasis() *
    quadrature.spaceDim() * quadrature.numBasis();
  for (int i=0; i < size; ++i)
    integrator._cellMatrix[i] = 1.23 + 1.2*i;
  integrator._resetCellMatrix();
  for (int i=0; i < size; ++i)
    CPPUNIT_ASSERT_EQUAL(0.0, integrator._cellMatrix[i]);
} // testResetCellMatrix

// ----------------------------------------------------------------------
// Set quadrature information.
void
pylith::feassemble::TestIntegrator::_initQuadrature(Quadrature1D* quadrature)
{ // _initQuadrature
  CPPUNIT_ASSERT(0 != quadrature);

  const int cellDim = 1;
  const int numBasis = 2;
  const int numQuadPts = 1;
  const int spaceDim = 1;
  const double vertices[] = { -1.0, 1.0 };
  const double basis[] = { 0.5, 0.5 };
  const double basisDeriv[] = { -0.5, 0.5 };
  const double quadPtsRef[] = { 0.0 };
  const double quadWts[] = { 2.0 };
  const double minJacobian = 1.0;

  quadrature->initialize(vertices, basis, basisDeriv, quadPtsRef, quadWts,
			 cellDim, numBasis, numQuadPts, spaceDim);
} // _initQuadrature


// End of file 
