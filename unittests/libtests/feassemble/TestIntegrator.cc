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
#include "pylith/feassemble/Quadrature1D.hh" // USES Quadrature1D
#include "pylith/utils/constdefs.h" // USES MAXDOUBLE

#include "spatialdata/spatialdb/GravityField.hh" // USES GravityField
#include "spatialdata/geocoords/CSCart.hh" // USES CSCart

#include <petscmat.h>

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

  Quadrature1D quadrature;
  const double minJacobian = 4.0;
  quadrature.minJacobian(minJacobian);
  
  ElasticityExplicit integrator;
  integrator.quadrature(&quadrature);

  CPPUNIT_ASSERT_EQUAL(minJacobian, integrator._quadrature->minJacobian());
} // testQuadrature

// ----------------------------------------------------------------------
// Test gravityField().
void
pylith::feassemble::TestIntegrator::testGravityField(void)
{ // testGravityField
  // Test gravity field by testing value of gravity vector.
  const int spaceDim = 3;
  const double gravityE[] = { 0.0, 0.0, -9.80665 };

  ElasticityImplicit integrator;
  spatialdata::spatialdb::GravityField gravityField;
  integrator.gravityField(&gravityField);

  spatialdata::geocoords::CSCart cs;
  cs.setSpaceDim(spaceDim);

  integrator._gravityField->open();
  double gravity[spaceDim];
  const double coords[] = { 1.0, 2.0, 3.0 };
  const int err = integrator._gravityField->query(gravity, spaceDim,
						  coords, spaceDim, &cs);
  CPPUNIT_ASSERT_EQUAL(0, err);
  integrator._gravityField->close();

  const double tolerance = 1.0e-06;
  for (int i=0; i < spaceDim; ++i)
    CPPUNIT_ASSERT_DOUBLES_EQUAL(gravityE[i], gravity[i], tolerance);
} // testGravityField

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
  const double basis[] = { 0.5, 0.5 };
  const double basisDeriv[] = { -0.5, 0.5 };
  const double quadPtsRef[] = { 0.0 };
  const double quadWts[] = { 2.0 };
  const double minJacobian = 1.0;

  quadrature->initialize(basis, basisDeriv, quadPtsRef, quadWts,
			 cellDim, numBasis, numQuadPts, spaceDim);
} // _initQuadrature


// End of file 
