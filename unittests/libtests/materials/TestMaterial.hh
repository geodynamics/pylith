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
  CPPUNIT_TEST( testTimeStep );
  CPPUNIT_TEST( testNeedNewJacobian );
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

  /// Test timeStep()
  void testTimeStep(void);

  /// Test needNewJacobian()
  void testNeedNewJacobian(void);

  /// Test initialize()
  void testInitialize(void);

  // PUBLIC METHODS /////////////////////////////////////////////////////
public :

  // Methods used in testing children of this class.

  /// Setup testing data.
  virtual
  void setUp(void);

  /// Tear down testing data.
  virtual
  void tearDown(void);

  /// Test dbToProperties().
  void testDBToProperties(void);

  /// Test dbValues().
  void testDBValues(void);

  /// Test _numProperties.
  void testProperties(void);

  // PROTECTED MEMBERS //////////////////////////////////////////////////
protected :

  Material* _material; ///< Object for testing
  MaterialData* _data; ///< Data for testing

}; // class TestMaterial

#endif // pylith_materials_testmaterial_hh

// End of file 
