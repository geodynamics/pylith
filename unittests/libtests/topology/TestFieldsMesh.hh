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
 * @file unittests/libtests/topology/TestFieldsMesh.hh
 *
 * @brief C++ unit testing for Fields<Mesh,Field>.
 */

#if !defined(pylith_topology_testfieldsmesh_hh)
#define pylith_topology_testfieldsmesh_hh

// Include directives ---------------------------------------------------
#include <cppunit/extensions/HelperMacros.h>

#include "pylith/topology/topologyfwd.hh" // forward declarations

// Forward declarations -------------------------------------------------
/// Namespace for pylith package
namespace pylith {
  namespace topology {
    class TestFieldsMesh;
  } // topology
} // pylith

// TestField -------------------------------------------------------------
/// C++ unit testing for Field.
class pylith::topology::TestFieldsMesh : public CppUnit::TestFixture
{ // class TestField

  // CPPUNIT TEST SUITE /////////////////////////////////////////////////
  CPPUNIT_TEST_SUITE( TestFieldsMesh );

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

}; // class TestFieldsMesh

#endif // pylith_topology_testfieldsmesh_hh


// End of file 
