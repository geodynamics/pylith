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
 * @file unittests/libtests/materials/TestElasticPlaneStress.hh
 *
 * @brief C++ TestElasticPlaneStress object
 *
 * C++ unit testing for ElasticPlaneStress.
 */

#if !defined(pylith_materials_testelasticplanestress_hh)
#define pylith_materials_testelasticplanestress_hh

#include "TestElasticMaterial.hh"

/// Namespace for pylith package
namespace pylith {
  namespace materials {
    class ElasticPlaneStress;
    class TestElasticPlaneStress;
    class ElasticPlaneStressData;
  } // materials
} // pylith

/// C++ unit testing for ElasticPlaneStress
class pylith::materials::TestElasticPlaneStress : public TestElasticMaterial
{ // class TestElasticPlaneStress

  // CPPUNIT TEST SUITE /////////////////////////////////////////////////
  CPPUNIT_TEST_SUITE( TestElasticPlaneStress );
  CPPUNIT_TEST( testDBToParameters );
  CPPUNIT_TEST( testDBValues );
  CPPUNIT_TEST( testParameters );
  CPPUNIT_TEST( testCalcDensity );
  CPPUNIT_TEST( testCalcStress );
  CPPUNIT_TEST( testCalcElasticConsts );
  CPPUNIT_TEST( testUpdateState );
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

  /// Test calcStress()
  void testCalcStress(void);

  /// Test calcElasticConsts()
  void testCalcElasticConsts(void);

  /// Test updateState()
  void testUpdateState(void);

}; // class TestElasticPlaneStress

#endif // pylith_materials_testelasticplanestress_hh


// End of file 
