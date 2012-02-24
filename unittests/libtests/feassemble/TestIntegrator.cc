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
// Copyright (c) 2010-2012 University of California, Davis
//
// See COPYING for license information.
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
  topology::Mesh mesh;

  CPPUNIT_ASSERT_EQUAL(pylith::PYLITH_MAXDOUBLE, integrator.stableTimeStep(mesh));
} // testStableTimeStep

// ----------------------------------------------------------------------
// Test isJacobianSymmetric().
void
pylith::feassemble::TestIntegrator::testIsJacobianSymmetric(void)
{ // testIsJacobianSymmetric
  ElasticityExplicit integrator;

  CPPUNIT_ASSERT_EQUAL(true, integrator.isJacobianSymmetric());

  integrator._isJacobianSymmetric = false;
  CPPUNIT_ASSERT_EQUAL(false, integrator.isJacobianSymmetric());
} // testIsJacobianSymmetric

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
// Test _lumpCellMatrix()
void
pylith::feassemble::TestIntegrator::testLumpCellMatrix(void)
{ // testLumpCellMatrix
  Quadrature<topology::Mesh> quadrature;
  _initQuadrature(&quadrature);

  ElasticityExplicit integrator;
  integrator.quadrature(&quadrature);

  integrator._initCellMatrix();
  integrator._initCellVector();
  
  const size_t sizeM = 
    quadrature.spaceDim() * quadrature.numBasis() *
    quadrature.spaceDim() * quadrature.numBasis();
  CPPUNIT_ASSERT_EQUAL(sizeM, integrator._cellMatrix.size());
  for (size_t i=0; i < sizeM; ++i)
    integrator._cellMatrix[i] = 1.23 + 1.2*i;
  integrator._lumpCellMatrix();

  const double tolerance = 1.0e-6;
  const int numBasis = quadrature.numBasis();
  const int spaceDim = quadrature.spaceDim();
  for (int iBasis=0; iBasis < numBasis; ++iBasis)
    for (int iDim=0; iDim < spaceDim; ++iDim) {
      double value = 0;
      const int index = (iBasis*spaceDim+iDim)*numBasis*spaceDim;
      for (int jBasis=0; jBasis < numBasis; ++jBasis)
	value += 1.23 + 1.2*(index+jBasis*spaceDim+iDim);
      CPPUNIT_ASSERT_DOUBLES_EQUAL(value, integrator._cellVector[iBasis*spaceDim+iDim], tolerance);
    } // for
} // testLumpCellMatrix

// ----------------------------------------------------------------------
// Test splitField().
void
pylith::feassemble::TestIntegrator::testSplitField(void)
{ // testSplitField

  topology::Mesh mesh;
  topology::Field<topology::Mesh> field(mesh);
  
  ElasticityExplicit integrator;
  integrator.splitField(&field);
  // Expect nothing to happen
} // testSplitField

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
