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

/**
 * @file unittests/libtests/materials/TestMaxwellIsotropic3D.hh
 *
 * @brief C++ TestMaxwellIsotropic3D object
 *
 * C++ unit testing for MaxwellIsotropic3D.
 */

#if !defined(pylith_materials_testmaxwellisotropic3d_hh)
#define pylith_materials_testmaxwellisotropic3d_hh

#include "TestElasticMaterial.hh"

/// Namespace for pylith package
namespace pylith {
  namespace materials {
    class TestMaxwellIsotropic3D;
    class MaxwellIsotropic3DElasticData;
    class MaxwellIsotropic3DTimeDepData;
  } // materials
} // pylith

/// C++ unit testing for MaxwellIsotropic3D
class pylith::materials::TestMaxwellIsotropic3D : public TestElasticMaterial
{ // class TestMaxwellIsotropic3D

  // CPPUNIT TEST SUITE /////////////////////////////////////////////////
  CPPUNIT_TEST_SUITE( TestMaxwellIsotropic3D );

  CPPUNIT_TEST( testDBToProperties );
  CPPUNIT_TEST( testNonDimProperties );
  CPPUNIT_TEST( testDimProperties );
  CPPUNIT_TEST( testDBToStateVars );
  CPPUNIT_TEST( testNonDimStateVars );
  CPPUNIT_TEST( testDimStateVars );

  CPPUNIT_TEST( test_calcDensity );

  CPPUNIT_TEST( testCalcStressElastic );
  CPPUNIT_TEST( testCalcStressTimeDep );
  CPPUNIT_TEST( testCalcElasticConstsElastic );
  CPPUNIT_TEST( testCalcElasticConstsTimeDep );
  CPPUNIT_TEST( testUpdatePropertiesElastic );
  CPPUNIT_TEST( testUpdatePropertiesTimeDep );
  CPPUNIT_TEST( test_stableTimeStepImplicit );

  CPPUNIT_TEST( testTimeStep );
  CPPUNIT_TEST( testUseElasticBehavior );
  CPPUNIT_TEST( testHasStateVars );


  CPPUNIT_TEST_SUITE_END();

  // PUBLIC METHODS /////////////////////////////////////////////////////
public :

  /// Setup testing data.
  void setUp(void);

  /// Test timeStep()
  void testTimeStep(void);

  /// Test useElasticBehavior()
  void testUseElasticBehavior(void);

  /// Test usesUpdateProperties()
  void testUsesUpdateProperties(void);

  /// Test calcStressElastic()
  void testCalcStressElastic(void);

  /// Test calcStressTimeDep()
  void testCalcStressTimeDep(void);

  /// Test calcElasticConstsElastic()
  void testCalcElasticConstsElastic(void);

  /// Test calcElasticConstsTimeDep()
  void testCalcElasticConstsTimeDep(void);

  /// Test updatePropertiesElastic()
  void testUpdatePropertiesElastic(void);

  /// Test updatePropertiesTimeDep()
  void testUpdatePropertiesTimeDep(void);

  /// Test _stableTimeStepImplicit()
  void test_stableTimeStepImplicit(void);

}; // class TestMaxwellIsotropic3D

#endif // pylith_materials_testmaxwellisotropic3d_hh


// End of file 
