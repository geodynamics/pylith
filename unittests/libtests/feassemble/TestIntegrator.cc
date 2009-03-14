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
#include "pylith/feassemble/ElasticityImplicit.hh" // USES ElasticityImplicit
#include "pylith/topology/Mesh.hh" // USES Mesh
#include "pylith/feassemble/Quadrature.hh" // USES Quadrature
#include "pylith/utils/constdefs.h" // USES MAXDOUBLE

#include "spatialdata/spatialdb/GravityField.hh" // USES GravityField

// ----------------------------------------------------------------------
CPPUNIT_TEST_SUITE_REGISTRATION( pylith::feassemble::TestIntegrator );

// ----------------------------------------------------------------------
// Test timeStep().
void
pylith::feassemble::TestIntegrator::testTimeStep(void)
{ // testTimeStep
  ElasticityExplicit integrator;
  const double dt = 1.2;
  integrator.timeStep(dt);

  CPPUNIT_ASSERT_EQUAL(dt, integrator._dt);
} // testTimeStep

// ----------------------------------------------------------------------
// Test stabletimeStep().
void
pylith::feassemble::TestIntegrator::testStableTimeStep(void)
{ // testStableTimeStep
  ElasticityExplicit integrator;

  CPPUNIT_ASSERT_EQUAL(pylith::PYLITH_MAXDOUBLE, integrator.stableTimeStep());
} // testStableTimeStep

// ----------------------------------------------------------------------
// Test quadrature().
void
pylith::feassemble::TestIntegrator::testQuadrature(void)
{ // testQuadrature
  // Since quadrature is cloned, test setting quadrature by testing
  // value of minJacobian

  Quadrature<topology::Mesh> quadrature;
  const double minJacobian = 4.0;
  quadrature.minJacobian(minJacobian);
  
  ElasticityExplicit integrator;
  integrator.quadrature(&quadrature);

  CPPUNIT_ASSERT_EQUAL(minJacobian, integrator._quadrature->minJacobian());
} // testQuadrature

// ----------------------------------------------------------------------
// Test normalizer().
void
pylith::feassemble::TestIntegrator::testNormalizer(void)
{ // testNormalizer
  const double lengthScale = 2.0;

  spatialdata::units::Nondimensional normalizer;
  normalizer.lengthScale(lengthScale);
  
  ElasticityExplicit integrator;
  integrator.normalizer(normalizer);

  CPPUNIT_ASSERT_EQUAL(lengthScale, integrator._normalizer->lengthScale());
} // testNormalizer

// ----------------------------------------------------------------------
// Test gravityField().
void
pylith::feassemble::TestIntegrator::testGravityField(void)
{ // testGravityField
  ElasticityImplicit integrator;
  spatialdata::spatialdb::GravityField gravityField;

  CPPUNIT_ASSERT(0 == integrator._gravityField);

  integrator.gravityField(&gravityField);
  CPPUNIT_ASSERT(0 != integrator._gravityField);
} // testGravityField

// ----------------------------------------------------------------------
// Test _initCellVector()
void
pylith::feassemble::TestIntegrator::testInitCellVector(void)
{ // testInitCellVector
  Quadrature<topology::Mesh> quadrature;
  _initQuadrature(&quadrature);

  ElasticityExplicit integrator;
  integrator.quadrature(&quadrature);

  integrator._initCellVector();
  
  const size_t size = 
    quadrature.spaceDim() * quadrature.numBasis();
  CPPUNIT_ASSERT_EQUAL(size, integrator._cellVector.size());
  for (size_t i=0; i < size; ++i)
    CPPUNIT_ASSERT_EQUAL(0.0, integrator._cellVector[i]);
} // testInitCellVector

// ----------------------------------------------------------------------
// Test _resetCellVector()
void
pylith::feassemble::TestIntegrator::testResetCellVector(void)
{ // testResetCellVector
  Quadrature<topology::Mesh> quadrature;
  _initQuadrature(&quadrature);

  ElasticityExplicit integrator;
  integrator.quadrature(&quadrature);

  integrator._initCellVector();
  
  const size_t size = 
    quadrature.spaceDim() * quadrature.numBasis();
  CPPUNIT_ASSERT_EQUAL(size, integrator._cellVector.size());
  for (size_t i=0; i < size; ++i)
    integrator._cellVector[i] = 1.4+2*i;
  integrator._resetCellVector();
  for (size_t i=0; i < size; ++i)
    CPPUNIT_ASSERT_EQUAL(0.0, integrator._cellVector[i]);
} // testResetCellVector

// ----------------------------------------------------------------------
// Test _initCellMatrix()
void
pylith::feassemble::TestIntegrator::testInitCellMatrix(void)
{ // testInitCellMatrix
  Quadrature<topology::Mesh> quadrature;
  _initQuadrature(&quadrature);

  ElasticityExplicit integrator;
  integrator.quadrature(&quadrature);

  integrator._initCellMatrix();
  
  const size_t size = 
    quadrature.spaceDim() * quadrature.numBasis() *
    quadrature.spaceDim() * quadrature.numBasis();
  CPPUNIT_ASSERT_EQUAL(size, integrator._cellMatrix.size());
  for (size_t i=0; i < size; ++i)
    CPPUNIT_ASSERT_EQUAL(0.0, integrator._cellMatrix[i]);
} // testInitCellMatrix

// ----------------------------------------------------------------------
// Test _resetCellMatrix()
void
pylith::feassemble::TestIntegrator::testResetCellMatrix(void)
{ // testResetCellMatrix
  Quadrature<topology::Mesh> quadrature;
  _initQuadrature(&quadrature);

  ElasticityExplicit integrator;
  integrator.quadrature(&quadrature);

  integrator._initCellMatrix();
  
  const size_t size = 
    quadrature.spaceDim() * quadrature.numBasis() *
    quadrature.spaceDim() * quadrature.numBasis();
  CPPUNIT_ASSERT_EQUAL(size, integrator._cellMatrix.size());
  for (size_t i=0; i < size; ++i)
    integrator._cellMatrix[i] = 1.23 + 1.2*i;
  integrator._resetCellMatrix();
  for (size_t i=0; i < size; ++i)
    CPPUNIT_ASSERT_EQUAL(0.0, integrator._cellMatrix[i]);
} // testResetCellMatrix

// ----------------------------------------------------------------------
// Set quadrature information.
void
pylith::feassemble::TestIntegrator::_initQuadrature(
				  Quadrature<topology::Mesh>* quadrature)
{ // _initQuadrature
  CPPUNIT_ASSERT(0 != quadrature);

  const int cellDim = 1;
  const int numBasis = 2;
  const int numQuadPts = 1;
  const int spaceDim = 1;
  const double basis[] = { 0.5, 0.5 };
  const double basisDeriv[] = { -0.5, 0.5 };
  const double quadPtsRef[] = { 0.0 };
  const double quadWts[] = { 2.0 };
  const double minJacobian = 1.0;

  quadrature->initialize(basis, numQuadPts, numBasis,
			 basisDeriv, numQuadPts, numBasis, cellDim,
			 quadPtsRef, numQuadPts, cellDim,
			 quadWts, numQuadPts,
			 spaceDim);
} // _initQuadrature


// End of file 
