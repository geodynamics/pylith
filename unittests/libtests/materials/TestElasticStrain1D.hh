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

}; // class TestElasticStrain1D

#endif // pylith_materials_testelasticstrain1d_hh


// End of file 
