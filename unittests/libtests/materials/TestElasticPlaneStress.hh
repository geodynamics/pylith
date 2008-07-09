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

  CPPUNIT_TEST( testDBToProperties );
  CPPUNIT_TEST( testDBValues );
  CPPUNIT_TEST( testProperties );
  CPPUNIT_TEST( test_calcDensity );
  CPPUNIT_TEST( test_calcStress );
  CPPUNIT_TEST( test_calcElasticConsts );
  CPPUNIT_TEST( test_stableTimeStepImplicit );

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

}; // class TestElasticPlaneStress

#endif // pylith_materials_testelasticplanestress_hh


// End of file 
