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
 * @file unittests/libtests/materials/TestIsotropicLinearElasticityPlaneStrain.hh
 *
 * @brief C++ TestIsotropicLinearElasticityPlaneStrain object
 *
 * C++ unit testing for IsotropicLinearElasticityPlaneStrain.
 */

#if !defined(pylith_materials_testisotropiclinearelasticityplanestrain_hh)
#define pylith_materials_testisotropiclinearelasticityplanestrain_hh

#include "TestMaterialNew.hh" // ISA TestMaterialNew

/// Namespace for pylith package
namespace pylith {
  namespace materials {
    class TestIsotropicLinearElasticityPlaneStrain;

    class TestIsotropicLinearElasticityPlaneStrain_Data;
  } // materials
} // pylith

/// C++ unit testing for IsotropicLinearElasticityPlaneStrain
class pylith::materials::TestIsotropicLinearElasticityPlaneStrain : public TestMaterialNew
{ // class TestIsotropicLinearElasticityPlaneStrain

  // CPPUNIT TEST SUITE /////////////////////////////////////////////////
  CPPUNIT_TEST_SUB_SUITE( TestIsotropicLinearElasticityPlaneStrain, TestMaterialNew );

  // Tests specific to this materials parameters.
  CPPUNIT_TEST( testUseInertia );
  CPPUNIT_TEST( testUseBodyForce );
  CPPUNIT_TEST( testUseInitialState );

  // Tests that explicitly depend on how details of this material.
  CPPUNIT_TEST( test_auxFieldsSetup );
  CPPUNIT_TEST( testGetAuxField );
  CPPUNIT_TEST( testIsJacobianSymmetric );

  CPPUNIT_TEST_SUITE_END();

  // PUBLIC METHODS /////////////////////////////////////////////////////
public :

  /// Setup testing data.
  void setUp(void);

  /// Deallocate testing data.
  void tearDown(void);

  /// Test useInertia().
  void testUseInertia(void);

  /// Test useBodyForce().
  void testUseBodyForce(void);

  /// Test useInitialState().
  void testUseInitialState(void);

  /// Test _auxFieldsSetup().
  void test_auxFieldsSetup(void);

  /// Test getAuxField().
  void testGetAuxField(void);

  /// Test IsJacobianSymmetric().
  void testIsJacobianSymmetric(void);

  // PROTECTED METHODS //////////////////////////////////////////////////
public :

  /** Get material.
   *
   * @returns Pointer to material.
   */
  MaterialNew* _material(void);

  /** Get test data.
   *
   * @returns Pointer to test data.
   */
  TestMaterialNew_Data* _data(void);

  /** Setup and populate solution field.
   *
   * @param[out] field Solution field to setup and populate.
   * @param[in] dbFilename Filename for spatial database with values for field.
   * @param[in] isClone True if field is a clone (don't need full setup).
   */
  void _setupSolutionField(pylith::topology::Field* field,
			   const char* dbFilename,
			   const bool isClone =false);

  // PROTECTED MEMBERS //////////////////////////////////////////////////
protected :

  IsotropicLinearElasticityPlaneStrain* _mymaterial; ///< Object for testing.
  TestIsotropicLinearElasticityPlaneStrain_Data* _mydata; ///< Data for testing.

}; // class TestIsotropicLinearElasticityPlaneStrain

#endif // pylith_materials_testisotropiclinearelasticityplanestrain_hh


// End of file 
