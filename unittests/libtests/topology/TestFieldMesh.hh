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
 * @file unittests/libtests/topology/TestFieldMesh.hh
 *
 * @brief C++ unit testing for Field.
 */

#if !defined(pylith_topology_testfieldmesh_hh)
#define pylith_topology_testfieldmesh_hh

// Include directives ---------------------------------------------------
#include <cppunit/extensions/HelperMacros.h>

#include "pylith/topology/topologyfwd.hh" // forward declarations

// Forward declarations -------------------------------------------------
/// Namespace for pylith package
namespace pylith {
  namespace topology {
    class TestFieldMesh;
  } // topology
} // pylith

// TestFieldMesh -------------------------------------------------------------
/// C++ unit testing for Field.
class pylith::topology::TestFieldMesh : public CppUnit::TestFixture
{ // class TestFieldMesh

  // CPPUNIT TEST SUITE /////////////////////////////////////////////////
  CPPUNIT_TEST_SUITE( TestFieldMesh );

  CPPUNIT_TEST( testConstructor );
  CPPUNIT_TEST( testMesh );
  CPPUNIT_TEST( testLabel );
  CPPUNIT_TEST( testVectorFieldType );
  CPPUNIT_TEST( testScale );
  CPPUNIT_TEST( testAddDimensionOkay );
  CPPUNIT_TEST( testSpaceDim );
  CPPUNIT_TEST( testChartSize );
  CPPUNIT_TEST( testSectionSize );
  CPPUNIT_TEST( testHasSection );
  CPPUNIT_TEST( testNewSectionPoints );
  CPPUNIT_TEST( testNewSectionPointsArray );
  CPPUNIT_TEST( testNewSectionDomain );
  CPPUNIT_TEST( testNewSectionField );
  CPPUNIT_TEST( testCloneSection );
  CPPUNIT_TEST( testClear );
  CPPUNIT_TEST( testAllocate );
  CPPUNIT_TEST( testZero );
  CPPUNIT_TEST( testZeroAll );
  CPPUNIT_TEST( testComplete );
  CPPUNIT_TEST( testCopy );
  CPPUNIT_TEST( testCopySubfield );
  CPPUNIT_TEST( testOperatorAdd );
  CPPUNIT_TEST( testDimensionalize );
  CPPUNIT_TEST( testView );
  CPPUNIT_TEST( testCreateScatter );
  CPPUNIT_TEST( testCreateScatterWithBC );
  CPPUNIT_TEST( testVector );
  CPPUNIT_TEST( testScatterLocalToGlobal );
  CPPUNIT_TEST( testScatterGlobalToLocal );
  CPPUNIT_TEST( testSplitDefault );
  CPPUNIT_TEST( testCloneSectionSplit );

  CPPUNIT_TEST_SUITE_END();

  // PUBLIC METHODS /////////////////////////////////////////////////////
public :

  /// Test constructor.
  void testConstructor(void);

  /// Test section().
  void testSection(void);

  /// Test mesh().
  void testMesh(void);

  /// Test label().
  void testLabel(void);

  /// Test vectorFieldType().
  void testVectorFieldType(void);

  /// Test scale().
  void testScale(void);

  /// Test addDimensionOkay().
  void testAddDimensionOkay(void);

  /// Test spaceDim().
  void testSpaceDim(void);

  /// Test chartSize().
  void testChartSize(void);

  /// Test sectionSize().
  void testSectionSize(void);

  /// Test hasSection().
  void testHasSection(void);

  /// Test newSection(points).
  void testNewSectionPoints(void);

  /// Test newSection(int_array).
  void testNewSectionPointsArray(void);

  /// Test newSection(domain).
  void testNewSectionDomain(void);

  /// Test newSection(field).
  void testNewSectionField(void);

  /// Test cloneSection().
  void testCloneSection(void);

  /// Test clear().
  void testClear(void);

  /// Test allocate().
  void testAllocate(void);

  /// Test zero().
  void testZero(void);

  /// Test zeroAll().
  void testZeroAll(void);

  /// Test complete().
  void testComplete(void);

  /// Test copy().
  void testCopy(void);

  /// Test copySubfield().
  void testCopySubfield(void);

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

  /// Test scatterLocalToGlobal().
  void testScatterLocalToGlobal(void);

  /// Test scatterGlobalToLocal().
  void testScatterGlobalToLocal(void);

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

}; // class TestFieldMesh

#endif // pylith_topology_testfieldmesh_hh


// End of file 
