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
 * @file unittests/libtests/materials/TestElasticStress1D.hh
 *
 * @brief C++ TestElasticStress1D object
 *
 * C++ unit testing for ElasticStress1D.
 */

#if !defined(pylith_materials_testelasticstress1d_hh)
#define pylith_materials_testelasticstress1d_hh

#include "TestElasticMaterial.hh"

/// Namespace for pylith package
namespace pylith {
  namespace materials {
    class ElasticStress1D;
    class TestElasticStress1D;
    class ElasticStress1DData;
  } // materials
} // pylith

/// C++ unit testing for ElasticStress1D
class pylith::materials::TestElasticStress1D : public TestElasticMaterial
{ // class TestElasticStress1D

  // CPPUNIT TEST SUITE /////////////////////////////////////////////////
  CPPUNIT_TEST_SUITE( TestElasticStress1D );

  CPPUNIT_TEST( testUsesUpdateState );
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

  /// Test usesUpdateState()
  void testUsesUpdateState();

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

}; // class TestElasticStress1D

#endif // pylith_materials_testelasticstress1d_hh


// End of file 
