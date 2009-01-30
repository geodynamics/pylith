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
 * @file unittests/libtests/topology/TestFieldBase.hh
 *
 * @brief C++ unit testing for FieldBase.
 */

#if !defined(pylith_topology_testfieldbase_hh)
#define pylith_topology_testfieldbase_hh

// Include directives ---------------------------------------------------
#include <cppunit/extensions/HelperMacros.h>

// Forward declarations -------------------------------------------------
/// Namespace for pylith package
namespace pylith {
  namespace topology {
    class TestFieldBase;

    class FieldBase;
  } // topology
} // pylith

// TestFieldBase -------------------------------------------------------------
/// C++ unit testing for FieldBase.
class pylith::topology::TestFieldBase : public CppUnit::TestFixture
{ // class TestFieldBase

  // CPPUNIT TEST SUITE /////////////////////////////////////////////////
  CPPUNIT_TEST_SUITE( TestFieldBase );

  CPPUNIT_TEST( testConstructor );
  CPPUNIT_TEST( testName );
  CPPUNIT_TEST( testVectorFieldType );
  CPPUNIT_TEST( testScale );
  CPPUNIT_TEST( testAddDimensionOkay );

  CPPUNIT_TEST_SUITE_END();

  // PUBLIC METHODS /////////////////////////////////////////////////////
public :

  /// Test constructor.
  void testConstructor(void);

  /// Test name().
  void testName(void);

  /// Test vectorFieldType().
  void testVectorFieldType(void);

  /// Test scale().
  void testScale(void);

  /// Test addDimensionOkay().
  void testAddDimensionOkay(void);

}; // class TestFieldBase

#endif // pylith_topology_testfieldbase_hh


// End of file 
