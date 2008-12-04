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
 * @file unittests/libtests/topology/TestFieldUniform.hh
 *
 * @brief C++ unit testing for FieldUniform.
 */

#if !defined(pylith_topology_testfielduniform_hh)
#define pylith_topology_testfielduniform_hh

// Include directives ---------------------------------------------------
#include <cppunit/extensions/HelperMacros.h>

// Forward declarations -------------------------------------------------
/// Namespace for pylith package
namespace pylith {
  namespace topology {
    class TestFieldUniform;

    class FieldUniform;
    class Mesh;
  } // topology
} // pylith

// TestFieldUniform -------------------------------------------------------------
/// C++ unit testing for FieldUniform.
class pylith::topology::TestFieldUniform : public CppUnit::TestFixture
{ // class TestFieldUniform

  // CPPUNIT TEST SUITE /////////////////////////////////////////////////
  CPPUNIT_TEST_SUITE( TestFieldUniform );

  CPPUNIT_TEST( testConstructor );
  CPPUNIT_TEST( testCreateSectionPoints );
  CPPUNIT_TEST( testCreateSectionChart );

  CPPUNIT_TEST_SUITE_END();

  // PUBLIC METHODS /////////////////////////////////////////////////////
public :

  /// Test constructor.
  void testConstructor(void);

  /// Test createSection() with points.
  void testCreateSectionPoints(void);

  /// Test createSection() with chart.
  void testCreateSectionChart(void);

// PRIVATE METHODS /////////////////////////////////////////////////////
private :

  /** Build mesh.
   *
   * @param mesh Finite-element mesh.
   */
  void _buildMesh(Mesh* mesh);

}; // class TestFieldUniform

#endif // pylith_topology_testfielduniform_hh


// End of file 
