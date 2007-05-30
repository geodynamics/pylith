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
  CPPUNIT_TEST( testUseElasticBehavior );
  CPPUNIT_TEST( testCalcDensity );
  CPPUNIT_TEST( testCalcStress );
  CPPUNIT_TEST( testCalcDerivElastic );
  CPPUNIT_TEST_SUITE_END();

  // PUBLIC METHODS /////////////////////////////////////////////////////
public :

  /// Test useElasticBehavior()
  void testUseElasticBehavior(void);

  /// Test calcDensity()
  void testCalcDensity(void);

  /// Test calcStress()
  void testCalcStress(void);

  /// Test calcDerivElastic()
  void testCalcDerivElastic(void);

  // PROTECTED METHODS //////////////////////////////////////////////////
protected :

  /** Test _calcDensity()
   *
   * @param material Pointer to material
   * @param data Data for testing material
   */
  void _testCalcDensity(ElasticMaterial* material,
			const ElasticMaterialData& data) const;

  /** Test _calcStress()
   *
   * @param material Pointer to material
   * @param data Data for testing material
   */
  void _testCalcStress(ElasticMaterial* material,
		       const ElasticMaterialData& data) const;

  /** Test _calcElasticConsts()
   *
   * @param material Pointer to material
   * @param data Data for testing material
   */
  void _testCalcElasticConsts(ElasticMaterial* material,
			      const ElasticMaterialData& data) const;

}; // class TestElasticMaterial

#endif // pylith_materials_testelasticmaterial_hh

// End of file 
