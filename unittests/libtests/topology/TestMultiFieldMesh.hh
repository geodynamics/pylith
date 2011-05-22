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
 * @file unittests/libtests/topology/TestMultiFieldMesh.hh
 *
 * @brief C++ unit testing for Field.
 */

#if !defined(pylith_topology_testmultifieldmesh_hh)
#define pylith_topology_testmultifieldmesh_hh

// Include directives ---------------------------------------------------
#include <cppunit/extensions/HelperMacros.h>

#include "pylith/topology/topologyfwd.hh" // forward declarations

// Forward declarations -------------------------------------------------
/// Namespace for pylith package
namespace pylith {
  namespace topology {
    class TestMultiFieldMesh;
  } // topology
} // pylith

// TestMultiFieldMesh -------------------------------------------------------------
/// C++ unit testing for Field.
class pylith::topology::TestMultiFieldMesh : public CppUnit::TestFixture
{ // class TestMultiFieldMesh

  // CPPUNIT TEST SUITE /////////////////////////////////////////////////
  CPPUNIT_TEST_SUITE( TestMultiFieldMesh );

  CPPUNIT_TEST( testConstructor );
  CPPUNIT_TEST( testDeallocate );

  CPPUNIT_TEST( testMesh );
  CPPUNIT_TEST( testSection );
  CPPUNIT_TEST( testLabel );
  CPPUNIT_TEST( testChartSize );
  CPPUNIT_TEST( testSectionSize );

  CPPUNIT_TEST( testHasField );
  CPPUNIT_TEST( testAdd );
  CPPUNIT_TEST( testGet );

  CPPUNIT_TEST( testFieldIndex );
  CPPUNIT_TEST( testFieldStartIndex );
  CPPUNIT_TEST( testFieldFiberDim );
  CPPUNIT_TEST( testFiberDim );

  CPPUNIT_TEST( testClear );
  CPPUNIT_TEST( testAllocate );
  CPPUNIT_TEST( testCloneSection );
  CPPUNIT_TEST( testComplete );

  CPPUNIT_TEST( testZero );
  CPPUNIT_TEST( testZeroAll );
  CPPUNIT_TEST( testCopy );
  CPPUNIT_TEST( testOperatorAdd );
  CPPUNIT_TEST( testDimensionalize );
  CPPUNIT_TEST( testView );

  CPPUNIT_TEST( testCreateScatter );
  CPPUNIT_TEST( testCreateScatterWithBC );
  CPPUNIT_TEST( testVector );
  CPPUNIT_TEST( testScatterSectionToVector );
  CPPUNIT_TEST( testScatterVectorToSection );
  CPPUNIT_TEST( testSplitDefault );
  CPPUNIT_TEST( testCloneSectionSplit );

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

  /// Test label().
  void testLabel(void);

  /// Test chartSize().
  void testChartSize(void);

  /// Test sectionSize().
  void testSectionSize(void);

  /// Test hasField().
  void testHadField(void);

  /// Test add().
  void testAdd(void);

  /// Test get().
  void testGet(void);

  /// Test fieldIndex().
  void testFieldIndex(void);

  /// Test fieldStartIndex().
  void testFieldStartIndex(void);

  /// Test fieldFiberDim().
  void testFieldFiberDim(void);

  /// Test fiberDim().
  void testFiberDim(void);

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

  /// Test zeroAll().
  void testZeroAll(void);

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

  /// Test splitDefault().
  void testSplitDefault(void);

  /// Test cloneSection() with split field.
  void testCloneSectionSplit(void);

// PRIVATE METHODS /////////////////////////////////////////////////////
private :

  /** Build mesh.
   *
   * @param mesh Finite-element mesh.
   */
  static
  void _buildMesh(Mesh* mesh);

}; // class TestMultiFieldMesh

#endif // pylith_topology_testmultifieldmesh_hh


// End of file 
