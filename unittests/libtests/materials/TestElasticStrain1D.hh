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
 * @file unittests/libtests/materials/TestElasticStrain1D.hh
 *
 * @brief C++ TestElasticStrain1D object
 *
 * C++ unit testing for ElasticStrain1D.
 */

#if !defined(pylith_materials_testelasticstrain1d_hh)
#define pylith_materials_testelasticstrain1d_hh

#include "TestElasticMaterial.hh"

/// Namespace for pylith package
namespace pylith {
  namespace materials {
    class ElasticStrain1D;
    class TestElasticStrain1D;
    class ElasticStrain1DData;
  } // materials
} // pylith

/// C++ unit testing for ElasticStrain1D
class pylith::materials::TestElasticStrain1D : public TestElasticMaterial
{ // class TestElasticStrain1D

  // CPPUNIT TEST SUITE /////////////////////////////////////////////////
  CPPUNIT_TEST_SUITE( TestElasticStrain1D );

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

}; // class TestElasticStrain1D

#endif // pylith_materials_testelasticstrain1d_hh


// End of file 
