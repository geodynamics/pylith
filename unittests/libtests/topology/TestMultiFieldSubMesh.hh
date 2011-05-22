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
// Copyright (c) 2010 University of California, Davis
//
// See COPYING for license information.
//
// ----------------------------------------------------------------------
//

/**
 * @file unittests/libtests/topology/TestMultiFieldSubMesh.hh
 *
 * @brief C++ unit testing for Field.
 */

#if !defined(pylith_topology_testmultifieldsubmesh_hh)
#define pylith_topology_testmultifieldsubmesh_hh

// Include directives ---------------------------------------------------
#include <cppunit/extensions/HelperMacros.h>

#include "pylith/topology/topologyfwd.hh" // forward declarations

// Forward declarations -------------------------------------------------
/// Namespace for pylith package
namespace pylith {
  namespace topology {
    class TestMultiFieldSubMesh;
  } // topology
} // pylith

// TestMultiFieldSubMesh -----------------------------------------------------
/// C++ unit testing for Field.
class pylith::topology::TestMultiFieldSubMesh : public CppUnit::TestFixture
{ // class TestMultiFieldSubMesh

  // CPPUNIT TEST SUITE /////////////////////////////////////////////////
  CPPUNIT_TEST_SUITE( TestMultiFieldSubMesh );

  CPPUNIT_TEST( testConstructor );
  CPPUNIT_TEST( testDeallocate );

  CPPUNIT_TEST( testMesh );
  CPPUNIT_TEST( testSection );

  CPPUNIT_TEST( testClear );
  CPPUNIT_TEST( testAllocate );
  CPPUNIT_TEST( testCloneSection );
  CPPUNIT_TEST( testComplete );

  CPPUNIT_TEST( testZero );
  CPPUNIT_TEST( testCopy );
  CPPUNIT_TEST( testOperatorAdd );
  CPPUNIT_TEST( testDimensionalize );
  CPPUNIT_TEST( testView );

  CPPUNIT_TEST( testCreateScatter );
  CPPUNIT_TEST( testCreateScatterWithBC );
  CPPUNIT_TEST( testVector );
  CPPUNIT_TEST( testScatterSectionToVector );
  CPPUNIT_TEST( testScatterVectorToSection );

  CPPUNIT_TEST_SUITE_END();

  // PUBLIC METHODS /////////////////////////////////////////////////////
public :

  /// Test constructor.
  void testConstructor(void);

  /// Test deallocate().
  void testDeallocate(void);

  /// Test mesh().
  void testMesh(void);

  /// Test section().
  void testSection(void);

  /// Test clear().
  void testClear(void);

  /// Test allocate().
  void testAllocate(void);

  /// Test cloneSection().
  void testCloneSection(void);

  /// Test complete().
  void testComplete(void);

  /// Test zero().
  void testZero(void);

  /// Test copy().
  void testCopy(void);

  /// Test operator+=().
  void testOperatorAdd(void);

  /// Test dimensionalize().
  void testDimensionalize(void);

  /// Test view().
  void testView(void);

  /// Test createScatter().
  void testCreateScatter(void);

  /// Test createScatterWithBC().
  void testCreateScatterWithBC(void);

  /// Test vector().
  void testVector(void);

  /// Test scatterSectionToVector().
  void testScatterSectionToVector(void);

  /// Test scatterVectorToSection().
  void testScatterVectorToSection(void);

// PRIVATE METHODS /////////////////////////////////////////////////////
private :

  /** Build mesh.
   *
   * @param mesh Finite-element mesh.
   * @param submesh Boundary mesh.
   */
  static
  void _buildMesh(Mesh* mesh,
		  SubMesh* submesh);

}; // class TestMultiFieldSubMesh

#endif // pylith_topology_testmultifieldsubmesh_hh


// End of file 
