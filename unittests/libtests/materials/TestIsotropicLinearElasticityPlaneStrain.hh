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

#include <cppunit/extensions/HelperMacros.h>

#include "pylith/materials/materialsfwd.hh" // forward declarations
#include "pylith/topology/topologyfwd.hh" // forward declarations

#include "spatialdata/spatialdb/SpatialDB.hh" // HOLDSA SpatialDB

/// Namespace for pylith package
namespace pylith {
  namespace materials {
    class TestIsotropicLinearElasticityPlaneStrain;

    class IsotropicLinearElasticityPlaneStrainData;
  } // materials
} // pylith

/// C++ unit testing for IsotropicLinearElasticityPlaneStrain
class pylith::materials::TestIsotropicLinearElasticityPlaneStrain : public CppUnit::TestFixture
{ // class TestIsotropicLinearElasticityPlaneStrain

  // CPPUNIT TEST SUITE /////////////////////////////////////////////////
  CPPUNIT_TEST_SUITE( TestIsotropicLinearElasticityPlaneStrain );

  // IsotropicLinearElasticityPlaneStrain
  CPPUNIT_TEST( testUseInertia );
  CPPUNIT_TEST( testUseBodyForce );
  CPPUNIT_TEST( testUseInitialState );
  CPPUNIT_TEST( test_auxFieldsSetup );
  CPPUNIT_TEST( test_setFEKernelsRHSResidual );
  CPPUNIT_TEST( test_setFEKernelsRHSJacobian );
  CPPUNIT_TEST( test_setFEKernelsLHSResidual );
  CPPUNIT_TEST( test_setFEKernelsLHSJacobianImplicit );
  CPPUNIT_TEST( test_setFEKernelsLHSJacobianExplicit );

  // Move to TestIntegratorPointwise
  CPPUNIT_TEST( testHasAuxField );
  CPPUNIT_TEST( testGetAuxField );
  CPPUNIT_TEST( testAuxFieldsDiscretization );
  CPPUNIT_TEST( testAuxFieldsDB );
  CPPUNIT_TEST( testNormalizer );
  CPPUNIT_TEST( testIsJacobianSymmetric );

  CPPUNIT_TEST( testVerifyConfiguration );
  CPPUNIT_TEST( testCheckConstraints );


  // Move to TestMaterialNew
  CPPUNIT_TEST( testDimension );
  CPPUNIT_TEST( testId );
  CPPUNIT_TEST( testLabel );
  CPPUNIT_TEST( testInitialize );

  CPPUNIT_TEST( testComputeResidual );
  CPPUNIT_TEST( testComputeRHSJacobian );
  CPPUNIT_TEST( testComputeLHSJacobianImplicit );
  CPPUNIT_TEST( testComputeLHSJacobianExplicit );
  CPPUNIT_TEST( testUpdateStateVars );

  CPPUNIT_TEST_SUITE_END();

  // PUBLIC METHODS /////////////////////////////////////////////////////
public :

  /// Setup testing data.
  virtual
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

  /// Test _setFEKernelsRHSResidual().
  void test_setFEKernelsRHSResidual(void);

  /// Test _setFEKernelsRHSJacobian().
  void test_setFEKernelsRHSJacobian(void);

  /// Test _setFEKernelsLHSResidual().
  void test_setFEKernelsLHSResidual(void);

  /// Test _setFEKernelsLHSJacobianImplicit().
  void test_setFEKernelsLHSJacobianImplicit(void);

  /// Test _setFEKernelsLHSJacobianExplicit().
  void test_setFEKernelsLHSJacobianExplicit(void);

  // IntegratorPointwise

  /// Test hasAuxField().
  void testHasAuxField(void);

  /// Test getAuxField().
  void testGetAuxField(void);

  /// Test auxFieldsDiscretization().
  void testAuxFieldsDiscretization(void);

  /// Test auxFieldsDB().
  void testAuxFieldsDB(void);

  /// Test normalizer().
  void testNormalizer(void);

  /// Test IsJacobianSymmetric().
  void testIsJacobianSymmetric(void);

  /// Test verifyConfiguration().
  void testVerifyConfiguration(void);

  /// Test checkConstraints().
  void testCheckConstraints(void);

  // MaterialNew

  /// Test dimension().
  void testDimension(void);

  /// Test id().
  void testId(void);

  /// Test label().
  void testLabel(void);

  /// Test initialize().
  void testInitialize(void);

  /// Test computeRHSResidual(), computeLHSResidual().
  void testComputeResidual(void);

  /// Test computeRHSJacobian().
  void testComputeRHSJacobian(void);

  /// Test computeLHSJacobianImplicit().
  void testComputeLHSJacobianImplicit(void);

  /// Test computeLHSJacobianExplicit().
  void testComputeLHSJacobianExplicit(void);

  /// Test updateStateVars().
  void testUpdateStateVars(void);

  // PROTECTED METHODS //////////////////////////////////////////////////
public :

  /// Do minimal initilaization of test data.
  void _initializeMin(void);

  /// Do full initilaization of test data.
  void _initializeFull(void);

  /** Setup and populate solution field.
   *
   * @param[out] field Solution field to setup and populate.
   * @param[in] dbFilename Filename for spatial database with values for field.
   * @param[in] isClone True if field is a clone (don't need full setup).
   */
  void _setupSolutionField(pylith::topology::Field* field,
			   const char* dbFilename,
			   const bool isClone =false);

  /** Set field to zero on the boundary.
   *
   * @param[out] field Field in which to set boundary values to zero.
   */
  void _zeroBoundary(pylith::topology::Field* field);

  // PROTECTED MEMBERS //////////////////////////////////////////////////
protected :

  IsotropicLinearElasticityPlaneStrain* _material; ///< Object for testing.
  IsotropicLinearElasticityPlaneStrainData* _data; ///< Data for testing.

  // TestMaterialNew
  topology::Mesh* _mesh; ///< Finite-element mesh.
  topology::Field* _solution1; ///< Solution field.
  topology::Field* _solution2; ///< Solution field.
  topology::Field* _solution1Dot; ///< Time derivative of solution field.
  topology::Field* _solution2Dot; ///< Time derivative of solution field.
  spatialdata::spatialdb::SimpleDB* _auxDB; ///< Spatial database with data for auxiliary fields.

}; // class TestIsotropicLinearElasticityPlaneStrain

#endif // pylith_materials_testisotropiclinearelasticityplanestrain_hh


// End of file 
