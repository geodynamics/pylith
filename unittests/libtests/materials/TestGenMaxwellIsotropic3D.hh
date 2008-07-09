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
 * @file unittests/libtests/materials/TestGenMaxwellIsotropic3D.hh
 *
 * @brief C++ TestGenMaxwellIsotropic3D object
 *
 * C++ unit testing for MaxwellIsotropic3D.
 */

#if !defined(pylith_materials_testgenmaxwellisotropic3d_hh)
#define pylith_materials_testgenmaxwellisotropic3d_hh

#include "TestElasticMaterial.hh"

/// Namespace for pylith package
namespace pylith {
  namespace materials {
    class GenMaxwellIsotropic3D;
    class TestGenMaxwellIsotropic3D;
    class GenMaxwellIsotropic3DElasticData;
    class GenMaxwellIsotropic3DTimeDepData;
  } // materials
} // pylith

/// C++ unit testing for GenMaxwellIsotropic3D
class pylith::materials::TestGenMaxwellIsotropic3D : public TestElasticMaterial
{ // class TestGenMaxwellIsotropic3D

  // CPPUNIT TEST SUITE /////////////////////////////////////////////////
  CPPUNIT_TEST_SUITE( TestGenMaxwellIsotropic3D );

  CPPUNIT_TEST( testDBToProperties );
  CPPUNIT_TEST( testDBValues );
  CPPUNIT_TEST( testProperties );
  CPPUNIT_TEST( test_calcDensity );

  CPPUNIT_TEST( testTimeStep );
  CPPUNIT_TEST( testUseElasticBehavior );
  CPPUNIT_TEST( testCalcStressElastic );
  CPPUNIT_TEST( testCalcStressTimeDep );
  CPPUNIT_TEST( testCalcElasticConstsElastic );
  CPPUNIT_TEST( testCalcElasticConstsTimeDep );
  CPPUNIT_TEST( testUsesUpdateProperties );
  CPPUNIT_TEST( testUpdatePropertiesElastic );
  CPPUNIT_TEST( testUpdatePropertiesTimeDep );
  CPPUNIT_TEST( test_stableTimeStepImplicit );

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

}; // class TestGenMaxwellIsotropic3D

#endif // pylith_materials_testgenmaxwellisotropic3d_hh


// End of file 
