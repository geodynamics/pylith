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
    class ElasticIsotropic3D;
    class TestElasticIsotropic3D;
    class ElasticIsotropic3DData;
  } // materials
} // pylith

/// C++ unit testing for ElasticIsotropic3D
class pylith::materials::TestElasticIsotropic3D : public TestElasticMaterial
{ // class TestElasticIsotropic3D

  // CPPUNIT TEST SUITE /////////////////////////////////////////////////
  CPPUNIT_TEST_SUITE( TestElasticIsotropic3D );

  CPPUNIT_TEST( testDBToProperties );
  CPPUNIT_TEST( testDBValues );
  CPPUNIT_TEST( testProperties );
  CPPUNIT_TEST( test_calcDensity );
  CPPUNIT_TEST( test_calcStress );
  CPPUNIT_TEST( test_calcElasticConsts );

  CPPUNIT_TEST( testUsesUpdateProperties );
  CPPUNIT_TEST( testUpdateProperties );

  CPPUNIT_TEST_SUITE_END();

  // PUBLIC METHODS /////////////////////////////////////////////////////
public :

  /// Setup testing data.
  void setUp(void);

  /// Test usesUpdateProperties().
  void testUsesUpdateProperties(void);

  /// Test updateProperties()
  void testUpdateProperties(void);

}; // class TestElasticIsotropic3D

#endif // pylith_materials_testelasticisotropic3d_hh


// End of file 
