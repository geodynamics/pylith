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

#include "TestIsotropicLinearElasticityPlaneStrain.hh" // Implementation of class methods

//#include "data/IsotropicLinearElasticityPlaneStrainData.hh" // USES IsotropicLinearElasticityPlaneStrainData

#include "pylith/materials/IsotropicLinearElasticityPlaneStrain.hh" // USES IsotropicLinearElasticityPlaneStrain

// ----------------------------------------------------------------------
CPPUNIT_TEST_SUITE_REGISTRATION( pylith::materials::TestIsotropicLinearElasticityPlaneStrain );

// ----------------------------------------------------------------------
// Setup testing data.
void
pylith::materials::TestIsotropicLinearElasticityPlaneStrain::setUp(void)
{ // setUp
  _material = new IsotropicLinearElasticityPlaneStrain();
#if 0
  _data = new IsotropicLinearElasticityPlaneStrainData();
  setupNormalizer();
#endif
} // setUp


// ----------------------------------------------------------------------
// Test useInertia().
void
pylith::materials::TestIsotropicLinearElasticityPlaneStrain::testUseInertia(void)
{ // testUseInertia
  CPPUNIT_ASSERT(false);
} // testUseInertia


// ----------------------------------------------------------------------
// Test useBodyForce().
void
pylith::materials::TestIsotropicLinearElasticityPlaneStrain::testUseBodyForce(void)
{ // testUseBodyForce
  CPPUNIT_ASSERT(false);
} // testUseBodyForce


// ----------------------------------------------------------------------
// Test preinitialize().
void
pylith::materials::TestIsotropicLinearElasticityPlaneStrain::testPreinitialize(void)
{ // testPreinitialize
  CPPUNIT_ASSERT(false);
} // testPreinitialize


// ----------------------------------------------------------------------
// Test _setFEKernels().
void
pylith::materials::TestIsotropicLinearElasticityPlaneStrain::test_setFEKernels(void)
{ // test_setFEKernels
  CPPUNIT_ASSERT(false);
} // test_setFEKernels


// ----------------------------------------------------------------------
// Test _dbToAuxFields().
void
pylith::materials::TestIsotropicLinearElasticityPlaneStrain::test_dbToAuxFields(void)
{ // test_dbToAuxFields
  CPPUNIT_ASSERT(false);
} // test_dbToAuxFields


// IntegratorPointwise ========================================


// ----------------------------------------------------------------------
// Test auxFields().
void
pylith::materials::TestIsotropicLinearElasticityPlaneStrain::testAuxFields(void)
{ // testAuxFields
  CPPUNIT_ASSERT(false);
} // testAuxFields


// ----------------------------------------------------------------------
// Test hasAuxField().
void
pylith::materials::TestIsotropicLinearElasticityPlaneStrain::testHasAuxField(void)
{ // testHasAuxField
  CPPUNIT_ASSERT(false);
} // testHaxAuxField


// ----------------------------------------------------------------------
// Test getAuxField().
void
pylith::materials::TestIsotropicLinearElasticityPlaneStrain::testGetAuxField(void)
{ // testGetAuxField
  CPPUNIT_ASSERT(false);
} // testGetAuxField


// ----------------------------------------------------------------------
// Test normalizer().
void
pylith::materials::TestIsotropicLinearElasticityPlaneStrain::testNormalizer(void)
{ // testNormalizer
  CPPUNIT_ASSERT(false);
} // testNormalizer


// ----------------------------------------------------------------------
// Test IsJacobianSymmetric().
void
pylith::materials::TestIsotropicLinearElasticityPlaneStrain::testIsJacobianSymmetric(void)
{ // testIsJacobianSymmetric
  CPPUNIT_ASSERT(false);
} // testIsJacobianSymmetric


// ----------------------------------------------------------------------
// Test verifyConfiguration().
void
pylith::materials::TestIsotropicLinearElasticityPlaneStrain::testVerifyConfiguration(void)
{ // testVerifyConfiguration
  CPPUNIT_ASSERT(false);
} // testVerifyConfiguration


// ----------------------------------------------------------------------
// Test initialize().
void
pylith::materials::TestIsotropicLinearElasticityPlaneStrain::testInitialize(void)
{ // testInitialize
  CPPUNIT_ASSERT(false);
} // testInitialize


// ----------------------------------------------------------------------
// Test integrateResidual().
void
pylith::materials::TestIsotropicLinearElasticityPlaneStrain::testIntegrateResidual(void)
{ // testIntegrateResidual
  CPPUNIT_ASSERT(false);
} // testIntegrateResidual


// ----------------------------------------------------------------------
// Test integrateJacobian().
void
pylith::materials::TestIsotropicLinearElasticityPlaneStrain::testIntegrateJacobian(void)
{ // testIntegrateJacobian
  CPPUNIT_ASSERT(false);
} // testIntegrateJacobian

// MaterialNew ========================================

// ----------------------------------------------------------------------
// Test dimension().
void
pylith::materials::TestIsotropicLinearElasticityPlaneStrain::testDimension(void)
{ // testDimension
  CPPUNIT_ASSERT(false);
} // testDimension


// ----------------------------------------------------------------------
// Test id().
void
pylith::materials::TestIsotropicLinearElasticityPlaneStrain::testId(void)
{ // testId
  CPPUNIT_ASSERT(false);
} // testId


// ----------------------------------------------------------------------
// Test label().
void
pylith::materials::TestIsotropicLinearElasticityPlaneStrain::testLabel(void)
{ // testLabel
  CPPUNIT_ASSERT(false);
} // testLabel


// ----------------------------------------------------------------------
// Test dbAuxFields().
void
pylith::materials::TestIsotropicLinearElasticityPlaneStrain::testDBAuxFields(void)
{ // testDBAuxFields
  CPPUNIT_ASSERT(false);
} // testDBAuxFields


// ----------------------------------------------------------------------
// Test _initializeAuxFieldsFromDB().
void
pylith::materials::TestIsotropicLinearElasticityPlaneStrain::test_initializeAuxFieldsFromDB(void)
{ // test_initializeAuxFieldsFromDB
  CPPUNIT_ASSERT(false);
} // test_initializeAuxFieldsFromDB


// ----------------------------------------------------------------------
// Test _nondimAuxFields().
void
pylith::materials::TestIsotropicLinearElasticityPlaneStrain::test_nondimAuxFields(void)
{ // test_nondimAuxFields
  CPPUNIT_ASSERT(false);
} // test_nondimAuxFields


// ----------------------------------------------------------------------
// Test _dimAuxFields().
void
pylith::materials::TestIsotropicLinearElasticityPlaneStrain::test_dimAuxFields(void)
{ // test_dimAuxFields
  CPPUNIT_ASSERT(false);
} // test_dimAuxFields


// End of file 
