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
    class MaxwellIsotropic3D;
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
  CPPUNIT_TEST( testDBToParameters );
  CPPUNIT_TEST( testDBValues );
  CPPUNIT_TEST( testParameters );
  CPPUNIT_TEST( testCalcDensity );
  CPPUNIT_TEST( testCalcStressElastic );
  CPPUNIT_TEST( testCalcStressTimeDep );
  CPPUNIT_TEST( testCalcElasticConstsElastic );
  CPPUNIT_TEST( testCalcElasticConstsTimeDep );
  CPPUNIT_TEST( testUpdateStateElastic );
  CPPUNIT_TEST( testUpdateStateTimeDep );
  CPPUNIT_TEST( testTimeStep );
  CPPUNIT_TEST_SUITE_END();

  // PUBLIC METHODS /////////////////////////////////////////////////////
public :

  /// Test DBValues()
  void testDBValues(void);

  /// Test parameters()
  void testParameters(void);

  /// Test _dbToParameters()
  void testDBToParameters(void);

  /// Test calcDensity()
  void testCalcDensity(void);

  /// Test calcStressElastic()
  void testCalcStressElastic(void);

  /// Test calcStressTimeDep()
  void testCalcStressTimeDep(void);

  /// Test calcElasticConstsElastic()
  void testCalcElasticConstsElastic(void);

  /// Test calcElasticConstsTimeDep()
  void testCalcElasticConstsTimeDep(void);

  /// Test updateStateElastic()
  void testUpdateStateElastic(void);

  /// Test updateStateTimeDep()
  void testUpdateStateTimeDep(void);

  /// Test timeStep()
  void testTimeStep(void);

}; // class TestMaxwellIsotropic3D

#endif // pylith_materials_testmaxwellisotropic3d_hh


// End of file 
