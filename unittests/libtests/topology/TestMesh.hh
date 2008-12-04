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
 * @file unittests/libtests/topology/TestMesh.hh
 *
 * @brief C++ unit testing for Mesh.
 */

#if !defined(pylith_topology_testmesh_hh)
#define pylith_topology_testmesh_hh

// Include directives ---------------------------------------------------
#include <cppunit/extensions/HelperMacros.h>

// Forward declarations -------------------------------------------------
/// Namespace for pylith package
namespace pylith {
  namespace topology {
    class TestMesh;

    class Mesh;
  } // topology
} // pylith

// TestMesh -------------------------------------------------------------
/// C++ unit testing for Mesh.
class pylith::topology::TestMesh : public CppUnit::TestFixture
{ // class TestMesh

  // CPPUNIT TEST SUITE /////////////////////////////////////////////////
  CPPUNIT_TEST_SUITE( TestMesh );

  CPPUNIT_TEST( testConstructor );
  CPPUNIT_TEST( testSieveMesh );
  CPPUNIT_TEST( testCoordsys );
  CPPUNIT_TEST( testInitialize );

  CPPUNIT_TEST_SUITE_END();

  // PUBLIC METHODS /////////////////////////////////////////////////////
public :

  /// Test constructor.
  void testConstructor(void);

  /// Test sieveMesh().
  void testSieveMesh(void);

  /// Test coordsys().
  void testCoordsys(void);

  /// Test initialize().
  void testInitialize(void);

}; // class TestMesh

#endif // pylith_topology_testmesh_hh


// End of file 
