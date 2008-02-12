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
 * @file unittests/libtests/materials/TestElasticMaterial.hh
 *
 * @brief C++ TestElasticMaterial object
 *
 * C++ unit testing for ElasticMaterial.
 */

#if !defined(pylith_materials_testelasticmaterial_hh)
#define pylith_materials_testelasticmaterial_hh

#include "TestMaterial.hh"

/// Namespace for pylith package
namespace pylith {
  namespace materials {
    class ElasticMaterial;
    class TestElasticMaterial;
    class ElasticMaterialData;
  } // materials
} // pylith

/// C++ unit testing for ElasticMaterial
class pylith::materials::TestElasticMaterial : public TestMaterial
{ // class TestElasticMaterial

  // CPPUNIT TEST SUITE /////////////////////////////////////////////////
  CPPUNIT_TEST_SUITE( TestElasticMaterial );

  CPPUNIT_TEST( testCalcDensity );
  CPPUNIT_TEST( testCalcStress );
  CPPUNIT_TEST( testCalcDerivElastic );

  CPPUNIT_TEST_SUITE_END();

  // PUBLIC METHODS /////////////////////////////////////////////////////
public :

  /// Test calcDensity()
  void testCalcDensity(void);

  /// Test calcStress()
  void testCalcStress(void);

  /// Test calcDerivElastic()
  void testCalcDerivElastic(void);

  // PUBLIC METHODS /////////////////////////////////////////////////////
public :

  // Methods used in testing children of this class.

  /// Setup testing data.
  virtual
  void setUp(void);

  /// Tear down testing data.
  virtual
  void tearDown(void);

  /// Test _calcDensity().
  void test_calcDensity(void);

  /// Test _calcStress().
  void test_calcStress(void);

  /// Test _calcElasticConsts().
  void test_calcElasticConsts(void);


  // PROTECTED MEMBERS //////////////////////////////////////////////////
protected :

  ElasticMaterial* _matElastic; ///< Test subject.
  ElasticMaterialData* _dataElastic; ///< Data for tests.

}; // class TestElasticMaterial

#endif // pylith_materials_testelasticmaterial_hh

// End of file 
