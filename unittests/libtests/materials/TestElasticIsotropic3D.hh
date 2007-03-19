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
  CPPUNIT_TEST( testDBToParameters );
  CPPUNIT_TEST( testDBValues );
  CPPUNIT_TEST( testParameters );
  CPPUNIT_TEST( testCalcDensity );
  CPPUNIT_TEST( testCalcElasticConsts );
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

  /// Test calcElasticConsts()
  void testCalcElasticConsts(void);

}; // class TestElasticIsotropic3D

#endif // pylith_materials_testelasticisotropic3d_hh


// End of file 
