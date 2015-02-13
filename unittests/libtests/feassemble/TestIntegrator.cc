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

#include "TestIntegrator.hh" // Implementation of class methods

#include "pylith/feassemble/ElasticityExplicit.hh" // USES ElasticityExplicit
#include "pylith/feassemble/ElasticityImplicit.hh" // USES ElasticityImplicit
#include "pylith/bc/Neumann.hh" // USES Neumann
#include "pylith/topology/Mesh.hh" // USES Mesh
#include "pylith/topology/Field.hh" // USES Field
#include "pylith/feassemble/Quadrature.hh" // USES Quadrature
#include "pylith/utils/constdefs.h" // USES MAXSCALAR

#include "spatialdata/spatialdb/GravityField.hh" // USES GravityField
#include "spatialdata/units/Nondimensional.hh" // USES Nondimensional

// ----------------------------------------------------------------------
CPPUNIT_TEST_SUITE_REGISTRATION( pylith::feassemble::TestIntegrator );

// ----------------------------------------------------------------------
// Test timeStep().
void
pylith::feassemble::TestIntegrator::testTimeStep(void)
{ // testTimeStep
  PYLITH_METHOD_BEGIN;

  ElasticityExplicit integrator;
  const PylithScalar dt = 1.2;
  integrator.timeStep(dt);

  CPPUNIT_ASSERT_EQUAL(dt, integrator._dt);

  PYLITH_METHOD_END;
} // testTimeStep

// ----------------------------------------------------------------------
// Test stabletimeStep().
void
pylith::feassemble::TestIntegrator::testStableTimeStep(void)
{ // testStableTimeStep
  PYLITH_METHOD_BEGIN;

  bc::Neumann integrator;
  topology::Mesh mesh;

  CPPUNIT_ASSERT_EQUAL(pylith::PYLITH_MAXSCALAR, integrator.stableTimeStep(mesh));

  PYLITH_METHOD_END;
} // testStableTimeStep

// ----------------------------------------------------------------------
// Test isJacobianSymmetric().
void
pylith::feassemble::TestIntegrator::testIsJacobianSymmetric(void)
{ // testIsJacobianSymmetric
  PYLITH_METHOD_BEGIN;

  ElasticityExplicit integrator;

  CPPUNIT_ASSERT_EQUAL(true, integrator.isJacobianSymmetric());

  integrator._isJacobianSymmetric = false;
  CPPUNIT_ASSERT_EQUAL(false, integrator.isJacobianSymmetric());

  PYLITH_METHOD_END;
} // testIsJacobianSymmetric

// ----------------------------------------------------------------------
// Test quadrature().
void
pylith::feassemble::TestIntegrator::testQuadrature(void)
{ // testQuadrature
  PYLITH_METHOD_BEGIN;

  // Since quadrature is cloned, test setting quadrature by testing
  // value of minJacobian

  Quadrature quadrature;
  const PylithScalar minJacobian = 4.0;
  quadrature.minJacobian(minJacobian);
  
  ElasticityExplicit integrator;
  integrator.quadrature(&quadrature);

  CPPUNIT_ASSERT_EQUAL(minJacobian, integrator._quadrature->minJacobian());

  PYLITH_METHOD_END;
} // testQuadrature

// ----------------------------------------------------------------------
// Test normalizer().
void
pylith::feassemble::TestIntegrator::testNormalizer(void)
{ // testNormalizer
  PYLITH_METHOD_BEGIN;

  const double lengthScale = 2.0;

  spatialdata::units::Nondimensional normalizer;
  normalizer.lengthScale(lengthScale);
  
  ElasticityExplicit integrator;
  integrator.normalizer(normalizer);

  CPPUNIT_ASSERT_EQUAL(lengthScale, integrator._normalizer->lengthScale());

  PYLITH_METHOD_END;
} // testNormalizer

// ----------------------------------------------------------------------
// Test gravityField().
void
pylith::feassemble::TestIntegrator::testGravityField(void)
{ // testGravityField
  PYLITH_METHOD_BEGIN;

  ElasticityImplicit integrator;
  spatialdata::spatialdb::GravityField gravityField;

  CPPUNIT_ASSERT(0 == integrator._gravityField);

  integrator.gravityField(&gravityField);
  CPPUNIT_ASSERT(integrator._gravityField);

  PYLITH_METHOD_END;
} // testGravityField

// ----------------------------------------------------------------------
// Test _initCellVector()
void
pylith::feassemble::TestIntegrator::testInitCellVector(void)
{ // testInitCellVector
  PYLITH_METHOD_BEGIN;

  Quadrature quadrature;
  _initQuadrature(&quadrature);

  ElasticityExplicit integrator;
  integrator.quadrature(&quadrature);

  integrator._initCellVector();
  
  const size_t size = 
    quadrature.spaceDim() * quadrature.numBasis();
  CPPUNIT_ASSERT_EQUAL(size, integrator._cellVector.size());
  for (size_t i=0; i < size; ++i)
    CPPUNIT_ASSERT_EQUAL(PylithScalar(0.0), integrator._cellVector[i]);

  PYLITH_METHOD_END;
} // testInitCellVector

// ----------------------------------------------------------------------
// Test _resetCellVector()
void
pylith::feassemble::TestIntegrator::testResetCellVector(void)
{ // testResetCellVector
  PYLITH_METHOD_BEGIN;

  Quadrature quadrature;
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
    CPPUNIT_ASSERT_EQUAL(PylithScalar(0.0), integrator._cellVector[i]);

  PYLITH_METHOD_END;
} // testResetCellVector

// ----------------------------------------------------------------------
// Test _initCellMatrix()
void
pylith::feassemble::TestIntegrator::testInitCellMatrix(void)
{ // testInitCellMatrix
  PYLITH_METHOD_BEGIN;

  Quadrature quadrature;
  _initQuadrature(&quadrature);

  ElasticityExplicit integrator;
  integrator.quadrature(&quadrature);

  integrator._initCellMatrix();
  
  const size_t size = 
    quadrature.spaceDim() * quadrature.numBasis() *
    quadrature.spaceDim() * quadrature.numBasis();
  CPPUNIT_ASSERT_EQUAL(size, integrator._cellMatrix.size());
  for (size_t i=0; i < size; ++i)
    CPPUNIT_ASSERT_EQUAL(PylithScalar(0.0), integrator._cellMatrix[i]);

  PYLITH_METHOD_END;
} // testInitCellMatrix

// ----------------------------------------------------------------------
// Test _resetCellMatrix()
void
pylith::feassemble::TestIntegrator::testResetCellMatrix(void)
{ // testResetCellMatrix
  PYLITH_METHOD_BEGIN;

  Quadrature quadrature;
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
    CPPUNIT_ASSERT_EQUAL(PylithScalar(0.0), integrator._cellMatrix[i]);

  PYLITH_METHOD_END;
} // testResetCellMatrix

// ----------------------------------------------------------------------
// Test _lumpCellMatrix()
void
pylith::feassemble::TestIntegrator::testLumpCellMatrix(void)
{ // testLumpCellMatrix
  PYLITH_METHOD_BEGIN;

  Quadrature quadrature;
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

  const PylithScalar tolerance = 1.0e-6;
  const int numBasis = quadrature.numBasis();
  const int spaceDim = quadrature.spaceDim();
  for (int iBasis=0; iBasis < numBasis; ++iBasis)
    for (int iDim=0; iDim < spaceDim; ++iDim) {
      PylithScalar value = 0;
      const int index = (iBasis*spaceDim+iDim)*numBasis*spaceDim;
      for (int jBasis=0; jBasis < numBasis; ++jBasis)
	value += 1.23 + 1.2*(index+jBasis*spaceDim+iDim);
      CPPUNIT_ASSERT_DOUBLES_EQUAL(value, integrator._cellVector[iBasis*spaceDim+iDim], tolerance);
    } // for

  PYLITH_METHOD_END;
} // testLumpCellMatrix

// ----------------------------------------------------------------------
// Set quadrature information.
void
pylith::feassemble::TestIntegrator::_initQuadrature(Quadrature* quadrature)
{ // _initQuadrature
  PYLITH_METHOD_BEGIN;

  CPPUNIT_ASSERT(quadrature);

  const int cellDim = 1;
  const int numBasis = 2;
  const int numQuadPts = 1;
  const int spaceDim = 1;
  const PylithScalar basis[numQuadPts*numBasis] = { 0.5, 0.5 };
  const PylithScalar basisDeriv[numQuadPts*numBasis*cellDim] = { -0.5, 0.5 };
  const PylithScalar quadPtsRef[numQuadPts*cellDim] = { 0.0 };
  const PylithScalar quadWts[numQuadPts] = { 2.0 };
  const PylithScalar minJacobian = 1.0;

  quadrature->initialize(basis, numQuadPts, numBasis,
			 basisDeriv, numQuadPts, numBasis, cellDim,
			 quadPtsRef, numQuadPts, cellDim,
			 quadWts, numQuadPts,
			 spaceDim);

  PYLITH_METHOD_END;
} // _initQuadrature


// End of file 
