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
// Copyright (c) 2010-2012 University of California, Davis
//
// See COPYING for license information.
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
  CPPUNIT_TEST( testCreateDMMesh );
  CPPUNIT_TEST( testDMMesh );
  CPPUNIT_TEST( testCoordsys );
  CPPUNIT_TEST( testDebug );
  CPPUNIT_TEST( testDimension );
  CPPUNIT_TEST( testComm );
  CPPUNIT_TEST( testNondimensionalize );

  CPPUNIT_TEST_SUITE_END();

  // PUBLIC METHODS /////////////////////////////////////////////////////
public :

  /// Test constructor.
  void testConstructor(void);

  /// Test createDMMesh().
  void testCreateDMMesh(void);

  /// Test dmMesh().
  void testDMMesh(void);

  /// Test coordsys().
  void testCoordsys(void);

  /// Test debug().
  void testDebug(void);

  /// Test dimension().
  void testDimension(void);

  /// Test comm().
  void testComm(void);

  /// Test nondimensionalize().
  void testNondimensionalize(void);

}; // class TestMesh

#endif // pylith_topology_testmesh_hh


// End of file 
