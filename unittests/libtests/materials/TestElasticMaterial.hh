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
    class TestElasticMaterial;
    class ElasticMaterialData;
  } // materials
} // pylith

/// C++ unit testing for ElasticMaterial
class pylith::materials::TestElasticMaterial : public TestMaterial
{ // class TestElasticMaterial

  // CPPUNIT TEST SUITE /////////////////////////////////////////////////
  CPPUNIT_TEST_SUITE( TestElasticMaterial );

  CPPUNIT_TEST( testDBInitialStress );
  CPPUNIT_TEST( testDBInitialStrain );
  CPPUNIT_TEST( testInitialize );
  CPPUNIT_TETS( testRetrievePropsAndVars );
  CPPUNIT_TEST( testCalcDensity );
  CPPUNIT_TEST( testCalcStress );
  CPPUNIT_TEST( testCalcDerivElastic );
  CPPUNIT_TEST( testUpdateStateVars );
  CPPUNIT_TEST( testStableTimeStepImplicit );
  CPPUNIT_TEST( testUseElasticBehavior );

  CPPUNIT_TEST_SUITE_END();

  // PUBLIC METHODS /////////////////////////////////////////////////////
public :

  /// Test dbInitialStress().
  void testDBInitialStress(void);

  /// Test dbInitialStrain().
  void testDBInitialStrain(void);

  /// Test initialize().
  void testInitialize(void);

  /// Test retrievePropsAndVars().
  void testRetrievePropsAndVars(void);

  /// Test calcDensity()
  void testCalcDensity(void);

  /// Test calcStress()
  void testCalcStress(void);

  /// Test calcDerivElastic()
  void testCalcDerivElastic(void);

  /// Test updateStateVars().
  void testUpdateStateVars(void);

  /// Test stableTimeStepImplicit().
  void testStableTimeStepImplicit(void);

  /// Test useElasticBehavior().
  void testUseElasticBehavior(void);

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

  /// Test _updateStateVars().
  void test_updateStateVars(void);

  /// Test _stableTimeStepImplicit().
  void test_stableTimeStepImplicit(void);


  // PROTECTED MEMBERS //////////////////////////////////////////////////
protected :

  ElasticMaterial* _matElastic; ///< Test subject.
  ElasticMaterialData* _dataElastic; ///< Data for tests.

}; // class TestElasticMaterial

#endif // pylith_materials_testelasticmaterial_hh

// End of file 
