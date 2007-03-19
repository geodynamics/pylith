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
 * @file unittests/libtests/materials/TestMaterial.hh
 *
 * @brief C++ TestMaterial object
 *
 * C++ unit testing for Material.
 */

#if !defined(pylith_materials_testmaterial_hh)
#define pylith_materials_testmaterial_hh

#include <cppunit/extensions/HelperMacros.h>

/// Namespace for pylith package
namespace pylith {
  namespace materials {
    class Material;
    class TestMaterial;
    class MaterialData;
  } // materials
} // pylith

/// C++ unit testing for Material
class pylith::materials::TestMaterial : public CppUnit::TestFixture
{ // class TestMaterial

  // CPPUNIT TEST SUITE /////////////////////////////////////////////////
  CPPUNIT_TEST_SUITE( TestMaterial );
  CPPUNIT_TEST( testDB );
  CPPUNIT_TEST( testID );
  CPPUNIT_TEST( testLabel );
  CPPUNIT_TEST( testInitialize );
  CPPUNIT_TEST_SUITE_END();

  // PUBLIC METHODS /////////////////////////////////////////////////////
public :

  /// Test db()
  void testDB(void);

  /// Test id()
  void testID(void);

  /// Test label()
  void testLabel(void);

  /// Test initialize()
  void testInitialize(void);

  // PROTECTED METHODS //////////////////////////////////////////////////
protected :

  /** Test dbToParameters().
   *
   * @param material Pointer to material
   * @param data Data for testing material
   */
  void _testDBToParameters(Material* material,
			   const MaterialData& data) const;

  /** Test _numDBValues and _dbValues()
   *
   * @param material Pointer to material
   * @param data Data for testing material
   */
  void _testDBValues(Material* material,
		     const MaterialData& data) const;

  /** Test _numParameters() and _parameterNames()
   *
   * @param material Pointer to material
   * @param data Data for testing material
   */
  void _testParameters(Material* material,
		       const MaterialData& data) const;

}; // class TestMaterial

#endif // pylith_materials_testmaterial_hh

// End of file 
