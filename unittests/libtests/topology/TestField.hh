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
 * @file unittests/libtests/topology/TestField.hh
 *
 * @brief C++ unit testing for Field.
 */

#if !defined(pylith_topology_testfield_hh)
#define pylith_topology_testfield_hh

// Include directives ---------------------------------------------------
#include <cppunit/extensions/HelperMacros.h>

// Forward declarations -------------------------------------------------
/// Namespace for pylith package
namespace pylith {
  namespace topology {
    class TestField;

    class Field;
    class Mesh;
  } // topology
} // pylith

// TestField -------------------------------------------------------------
/// C++ unit testing for Field.
class pylith::topology::TestField : public CppUnit::TestFixture
{ // class TestField

  // CPPUNIT TEST SUITE /////////////////////////////////////////////////
  CPPUNIT_TEST_SUITE( TestField );

  CPPUNIT_TEST( testConstructor );
  CPPUNIT_TEST( testSection );
  CPPUNIT_TEST( testSpaceDim );
  CPPUNIT_TEST( testNewSection );
  CPPUNIT_TEST( testNewSectionPoints );
  CPPUNIT_TEST( testNewSectionDomain );
  CPPUNIT_TEST( testNewSectionChart );
  CPPUNIT_TEST( testNewSectionField );
  CPPUNIT_TEST( testClear );
  CPPUNIT_TEST( testAllocate );
  CPPUNIT_TEST( testZero );
  CPPUNIT_TEST( testComplete );
  CPPUNIT_TEST( testCopy );
  CPPUNIT_TEST( testOperatorAdd );
  CPPUNIT_TEST( testDimensionalize );
  CPPUNIT_TEST( testView );

  CPPUNIT_TEST_SUITE_END();

  // PUBLIC METHODS /////////////////////////////////////////////////////
public :

  /// Test constructor.
  void testConstructor(void);

  /// Test section().
  void testSection(void);

  /// Test spaceDim().
  void testSpaceDim(void);

  /// Test newSection().
  void testNewSection(void);

  /// Test newSection(points).
  void testNewSectionPoints(void);

  /// Test newSection(domain).
  void testNewSectionDomain(void);

  /// Test newSection(chart).
  void testNewSectionChart(void);

  /// Test newSection(field).
  void testNewSectionField(void);

  /// Test clear().
  void testClear(void);

  /// Test allocate().
  void testAllocate(void);

  /// Test zero().
  void testZero(void);

  /// Test complete().
  void testComplete(void);

  /// Test copy().
  void testCopy(void);

  /// Test operator+=().
  void testOperatorAdd(void);

  /// Test dimensionalize().
  void testDimensionalize(void);

  /// Test view().
  void testView(void);

// PRIVATE METHODS /////////////////////////////////////////////////////
private :

  /** Build mesh.
   *
   * @param mesh Finite-element mesh.
   */
  void _buildMesh(Mesh* mesh);

}; // class TestField

#endif // pylith_topology_testfield_hh


// End of file 
