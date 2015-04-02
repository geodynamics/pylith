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
// Copyright (c) 2010-2015 University of California, Davis
//
// See COPYING for license information.
//
// ----------------------------------------------------------------------
//

/**
 * @file unittests/libtests/topology/TestFieldsSubMesh.hh
 *
 * @brief C++ unit testing for Fields<Mesh,Field>.
 */

#if !defined(pylith_topology_testfieldssubmesh_hh)
#define pylith_topology_testfieldssubmesh_hh

// Include directives ---------------------------------------------------
#include <cppunit/extensions/HelperMacros.h>

// Forward declarations -------------------------------------------------
/// Namespace for pylith package
namespace pylith {
  namespace topology {
    class TestFieldsSubMesh;

    class Mesh;
  } // topology
} // pylith

// TestField -------------------------------------------------------------
/// C++ unit testing for Field.
class pylith::topology::TestFieldsSubMesh : public CppUnit::TestFixture
{ // class TestField

  // CPPUNIT TEST SUITE /////////////////////////////////////////////////
  CPPUNIT_TEST_SUITE( TestFieldsSubMesh );

  CPPUNIT_TEST( testConstructor );
  CPPUNIT_TEST( testAdd );
  CPPUNIT_TEST( testAddDomain );
  CPPUNIT_TEST( testDelete );
  CPPUNIT_TEST( testGet );
  CPPUNIT_TEST( testGetConst );
  CPPUNIT_TEST( testHasField );
  CPPUNIT_TEST( testCopyLayout );

  CPPUNIT_TEST_SUITE_END();

  // PUBLIC METHODS /////////////////////////////////////////////////////
public :

  /// Setup test case.
  void setUp(void);

  /// Tear down test case.
  void tearDown(void);

  /// Test constructor.
  void testConstructor(void);

  /// Test add().
  void testAdd(void);

  /// Test add(domain).
  void testAddDomain(void);

  /// Test delete().
  void testDelete(void);

  /// Test get().
  void testGet(void);

  /// Test get() for const Fields.
  void testGetConst(void);

  /// Test hasField().
  void testHasField(void);

  /// Test copyLayout(domain).
  void testCopyLayout(void);

// PRIVATE MEMBERS /////////////////////////////////////////////////////
private :

  Mesh* _mesh;
  Mesh* _submesh;

}; // class TestFieldsSubMesh

#endif // pylith_topology_testfieldssubmesh_hh


// End of file 
