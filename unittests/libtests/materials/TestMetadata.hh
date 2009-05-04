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
 * @file unittests/libtests/materials/TestMetadata.hh
 *
 * @brief C++ TestMetadata object
 *
 * C++ unit testing for Material.
 */

#if !defined(pylith_materials_testmetadata_hh)
#define pylith_materials_testmetadata_hh

#include <cppunit/extensions/HelperMacros.h>

#include "pylith/materials/materialsfwd.hh" // forward declarations

/// Namespace for pylith package
namespace pylith {
  namespace materials {
    class TestMetadata;
  } // materials
} // pylith

/// C++ unit testing for Material
class pylith::materials::TestMetadata : public CppUnit::TestFixture
{ // class TestMetadata

  // CPPUNIT TEST SUITE /////////////////////////////////////////////////
  CPPUNIT_TEST_SUITE( TestMetadata );

  CPPUNIT_TEST( testConstructor );
  CPPUNIT_TEST( testCopyConstructor );
  CPPUNIT_TEST( testProperties );
  CPPUNIT_TEST( testStateVars );
  CPPUNIT_TEST( testFiberDim );
  CPPUNIT_TEST( testFieldType );
  CPPUNIT_TEST( testDBProperties );
  CPPUNIT_TEST( testDBStateVars );

  CPPUNIT_TEST_SUITE_END();

  // PUBLIC METHODS /////////////////////////////////////////////////////
public :

  /// Setup test data.
  void setUp(void);

  /// Tear down test data.
  void tearDown(void);

  /// Test constructor.
  void testConstructor(void);

  /// Test copy constructor.
  void testCopyConstructor(void);

  /// Test properties().
  void testProperties(void);

  /// Test stateVars().
  void testStateVars(void);

  /// Test fiberDim().
  void testFiberDim(void);

  /// Test fieldType().
  void testFieldType(void);

  /// Test dbProperties() and numDBProperties().
  void testDBProperties(void);

  /// Test dbStateVars() and numDBStateVars().
  void testDBStateVars(void);

  // PRIVATE MEMBERS ////////////////////////////////////////////////////
private :

  Metadata* _metadata; ///< Object for testing

}; // class TestMetadata

#endif // pylith_materials_testmetadata_hh

// End of file 
