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
 * @file unittests/libtests/materials/TestElasticIsotropic3D.hh
 *
 * @brief C++ TestElasticIsotropic3D object
 *
 * C++ unit testing for ElasticIsotropic3D.
 */

#if !defined(pylith_materials_testelasticisotropic3d_hh)
#define pylith_materials_testelasticisotropic3d_hh

#include "TestElasticMaterial.hh"

/// Namespace for pylith package
namespace pylith {
  namespace materials {
    class TestElasticIsotropic3D;
  } // materials
} // pylith

/// C++ unit testing for ElasticIsotropic3D
class pylith::materials::TestElasticIsotropic3D : public TestElasticMaterial
{ // class TestElasticIsotropic3D

  // CPPUNIT TEST SUITE /////////////////////////////////////////////////
  CPPUNIT_TEST_SUITE( TestElasticIsotropic3D );

  CPPUNIT_TEST( testDimension );
  CPPUNIT_TEST( testTensorSize );
  CPPUNIT_TEST( testDBToProperties );
  CPPUNIT_TEST( testNonDimProperties );
  CPPUNIT_TEST( testDimProperties );
  CPPUNIT_TEST( testDBToStateVars );
  CPPUNIT_TEST( testNonDimStateVars );
  CPPUNIT_TEST( testDimStateVars );
  CPPUNIT_TEST( testStableTimeStepImplicit );
  CPPUNIT_TEST( test_calcDensity );
  CPPUNIT_TEST( test_calcStress );
  CPPUNIT_TEST( test_calcElasticConsts );
  CPPUNIT_TEST( test_updateStateVars );
  CPPUNIT_TEST( test_stableTimeStepImplicit );

  CPPUNIT_TEST_SUITE_END();

  // PUBLIC METHODS /////////////////////////////////////////////////////
public :

  /// Setup testing data.
  void setUp(void);

  /// Test stableTimeStepImplicit().
  void testStableTimeStepImplicit(void);

}; // class TestElasticIsotropic3D

#endif // pylith_materials_testelasticisotropic3d_hh


// End of file 
