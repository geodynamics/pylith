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

/**
 * @file unittests/libtests/materials/TestIostropicLinearElasticityPlaneStrain.hh
 *
 * @brief C++ TestIostropicLinearElasticityPlaneStrain object
 *
 * C++ unit testing for IostropicLinearElasticityPlaneStrain.
 */

#if !defined(pylith_materials_testisitropiclinearelasticityplanestrain_hh)
#define pylith_materials_testisotropiclinearelasticityplanestrain_hh

#include <cppunit/extensions/HelperMacros.h>

#include "pylith/materials/materialsfwd.hh" // forward declarations

/// Namespace for pylith package
namespace pylith {
  namespace materials {
    class TestIsotropicLinearElasticityPlaneStrain;
  } // materials
} // pylith

/// C++ unit testing for IsotropicLinearElasticityPlaneStrain
class pylith::materials::TestIsotropicLinearElasticityPlaneStrain : public CppUnit::TestFixture
{ // class TestIsotropicLinearElasticityPlaneStrain

  // CPPUNIT TEST SUITE /////////////////////////////////////////////////
  CPPUNIT_TEST_SUITE( TestIsotropicLinearElasticityPlaneStrain );

  CPPUNIT_TEST( testUseInertia );
  CPPUNIT_TEST( testUseBodyForce );

  CPPUNIT_TEST( testPreinitialize );
  CPPUNIT_TEST( test_setFEKernels );
  CPPUNIT_TEST( test_dbToAuxFields );

  // Move to TestIntegratorPointwise
  CPPUNIT_TEST( testAuxFields );
  CPPUNIT_TEST( testHasAuxField );
  CPPUNIT_TEST( testGetAuxField );
  CPPUNIT_TEST( testNormalizer );
  CPPUNIT_TEST( testIsJacobianSymmetric );

  CPPUNIT_TEST( testVerifyConfiguration );
  CPPUNIT_TEST( testInitialize );
  CPPUNIT_TEST( testIntegrateResidual );
  CPPUNIT_TEST( testIntegrateJacobian );

  // Move to TestMaterialNew
  CPPUNIT_TEST( testDimension );
  CPPUNIT_TEST( testId );
  CPPUNIT_TEST( testLabel );
  CPPUNIT_TEST( testDBAuxFields );

  CPPUNIT_TEST( test_initializeAuxFieldsFromDB );
  CPPUNIT_TEST( test_nondimAuxFields );
  CPPUNIT_TEST( test_dimAuxFields );

  CPPUNIT_TEST_SUITE_END();

  // PUBLIC METHODS /////////////////////////////////////////////////////
public :

  /// Setup testing data.
  void setUp(void);

  /// Test dimension().
  void testDimension(void);

  /// Test id().
  void testId(void);

  /// Test label().
  void testLabel(void);

  /// Test useInertia().
  void testUseInertia(void);

  /// Test useBodyForce().
  void testUseBodyForce(void);

  /// Test preinitialize().
  void testPreinitialize(void);

  /// Test _setFEKernels().
  void test_setFEKernels(void);

  /// Test _dbToAuxFields().
  void test_dbToAuxFields(void);

  // IntegratorPointwise

  /// Test auxFields().
  void testAuxFields(void);

  /// Test hasAuxField().
  void testHasAuxField(void);

  /// Test getAuxField().
  void testGetAuxField(void);

  /// Test normalizer().
  void testNormalizer(void);

  /// Test IsJacobianSymmetric().
  void testIsJacobianSymmetric(void);

  /// Test verifyConfiguration().
  void testVerifyConfiguration(void);

  /// Test initialize().
  void testInitialize(void);

  /// Test integrateResidual().
  void testIntegrateResidual(void);

  /// Test integrateJacobian().
  void testIntegrateJacobian(void);

  // MaterialNew

  /// Test dbAuxFields().
  void testDBAuxFields(void);

  /// Test _initializeAuxFieldsFromDB().
  void test_initializeAuxFieldsFromDB(void);

  /// Test _nondimAuxFields().
  void test_nondimAuxFields(void);

  /// Test _dimAuxFields().
  void test_dimAuxFields(void);

  // PROTECTED MEMBERS //////////////////////////////////////////////////
protected :

  IsotropicLinearElasticityPlaneStrain* _material; ///< Object for testing
  //MaterialNewData* _data; ///< Data for testing

}; // class TestIsotropicLinearElasticityPlaneStrain

#endif // pylith_materials_testisotropiclinearelasticityplanestrain_hh


// End of file 
