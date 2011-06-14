// -*- C++ -*-
//
// ----------------------------------------------------------------------
//
// Brad T. Aagaard, U.S. Geological Survey
// Charles A. Williams, GNS Science
// Matthew G. Knepley, University of Chicago
//
// This code was developed as part of the Computational Infrastructure
// for Geodynamics (http://geodynamics.org).
//
// Copyright (c) 2010-2011 University of California, Davis
//
// See COPYING for license information.
//
// ----------------------------------------------------------------------
//

/**
 * @file unittests/libtests/topology/TestSubMesh.hh
 *
 * @brief C++ unit testing for Mesh.
 */

#if !defined(pylith_topology_testsubmesh_hh)
#define pylith_topology_testsubmesh_hh

// Include directives ---------------------------------------------------
#include <cppunit/extensions/HelperMacros.h>

// Forward declarations -------------------------------------------------
/// Namespace for pylith package
namespace pylith {
  namespace topology {
    class TestSubMesh;

    class Mesh;
  } // topology
} // pylith

// TestSubMesh -------------------------------------------------------------
/// C++ unit testing for Mesh.
class pylith::topology::TestSubMesh : public CppUnit::TestFixture
{ // class TestSubMesh

  // CPPUNIT TEST SUITE /////////////////////////////////////////////////
  CPPUNIT_TEST_SUITE( TestSubMesh );

  CPPUNIT_TEST( testConstructor );
  CPPUNIT_TEST( testConstructorMesh );
  CPPUNIT_TEST( testSieveMesh );
  CPPUNIT_TEST( testCreateSubMesh );
  CPPUNIT_TEST( testCoordsys );
  CPPUNIT_TEST( testDebug );
  CPPUNIT_TEST( testDimension );
  CPPUNIT_TEST( testComm );
  CPPUNIT_TEST( testInitialize );

  CPPUNIT_TEST_SUITE_END();

  // PUBLIC METHODS /////////////////////////////////////////////////////
public :

  /// Test constructor.
  void testConstructor(void);

  /// Test constructor w/mesh.
  void testConstructorMesh(void);

  /// Test sieveMesh().
  void testSieveMesh(void);

  /// Test createSubMesh.
  void testCreateSubMesh(void);

  /// Test coordsys().
  void testCoordsys(void);

  /// Test debug().
  void testDebug(void);

  /// Test dimension().
  void testDimension(void);

  /// Test comm().
  void testComm(void);

  /// Test initialize().
  void testInitialize(void);

// PRIVATE METHODS /////////////////////////////////////////////////////
private :

  /** Build mesh.
   *
   * @param mesh Finite-element mesh.
   */
  static
  void _buildMesh(Mesh* mesh);

}; // class TestSubMesh

#endif // pylith_topology_testsubmesh_hh


// End of file 
