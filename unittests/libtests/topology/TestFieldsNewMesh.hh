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
// Copyright (c) 2010-2013 University of California, Davis
//
// See COPYING for license information.
//
// ----------------------------------------------------------------------
//

/**
 * @file unittests/libtests/topology/TestFieldsNewMesh.hh
 *
 * @brief C++ unit testing for FieldsNew<Mesh>.
 */

#if !defined(pylith_topology_testfieldsnewmesh_hh)
#define pylith_topology_testfieldsnewmesh_hh

// Include directives ---------------------------------------------------
#include <cppunit/extensions/HelperMacros.h>

#include "pylith/topology/topologyfwd.hh" // forward declarations

// Forward declarations -------------------------------------------------
/// Namespace for pylith package
namespace pylith {
  namespace topology {
    class TestFieldsNewMesh;
  } // topology
} // pylith

// TestField -------------------------------------------------------------
/// C++ unit testing for Field.
class pylith::topology::TestFieldsNewMesh : public CppUnit::TestFixture
{ // class TestField

  // CPPUNIT TEST SUITE /////////////////////////////////////////////////
  CPPUNIT_TEST_SUITE( TestFieldsNewMesh );

  CPPUNIT_TEST( testConstructor );
  CPPUNIT_TEST( testHasField );
  CPPUNIT_TEST( testAdd );
  CPPUNIT_TEST( testAllocateSequence );
  CPPUNIT_TEST( testAllocateArray );
  CPPUNIT_TEST( testAllocateDomain );
  CPPUNIT_TEST( testGet );
  CPPUNIT_TEST( testGetConst );
  CPPUNIT_TEST( testMesh );
  CPPUNIT_TEST( testFiberDim );
  CPPUNIT_TEST( testSectionIndex );
  CPPUNIT_TEST( testSectionFiberDim );
  CPPUNIT_TEST( testComplete );
  CPPUNIT_TEST( testFieldNames );
  CPPUNIT_TEST( testView );

  CPPUNIT_TEST_SUITE_END();

  // PUBLIC METHODS /////////////////////////////////////////////////////
public :

  /// Setup test case.
  void setUp(void);

  /// Tear down test case.
  void tearDown(void);

  /// Test constructor.
  void testConstructor(void);

  /// Test hasField().
  void testHasField(void);

  /// Test add().
  void testAdd(void);

  /// Test allocate(label_sequence).
  void testAllocateSequence(void);

  /// Test allocate(int_array).
  void testAllocateArray(void);

  /// Test allocate(domain).
  void testAllocateDomain(void);

  /// Test getField().
  void testGet(void);

  /// Test getField() for const FieldsNew.
  void testGetConst(void);

  /// Test mesh().
  void testMesh(void);

  /// Test fiberDim().
  void testFiberDim(void);

  /// Test sectionIndex().
  void testSectionIndex(void);

  /// Test sectionFiberDim().
  void testSectionFiberDim(void);

  /// Test complete().
  void testComplete(void);

  /// Test fieldNames().
  void testFieldNames(void);

  /// Test view().
  void testView(void);

// PRIVATE MEMBERS /////////////////////////////////////////////////////
private :

  Mesh* _mesh;

}; // class TestFieldsNewMesh

#endif // pylith_topology_testfieldsnewmesh_hh


// End of file 
