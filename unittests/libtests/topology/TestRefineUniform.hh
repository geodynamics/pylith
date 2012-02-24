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
 * @file unittests/libtests/topology/TestRefineUniform.hh
 *
 * @brief C++ TestRefineUniform object
 *
 * C++ unit testing for RefineUniform.
 */

#if !defined(pylith_topology_testrefineuniform_hh)
#define pylith_topology_testrefineuniform_hh

// Include directives ---------------------------------------------------
#include <cppunit/extensions/HelperMacros.h>

#include "pylith/topology/topologyfwd.hh" // USES Mesh

// Forward declarations -------------------------------------------------
/// Namespace for pylith package
namespace pylith {
  namespace topology {
    class TestRefineUniform;

    class MeshDataCohesive; // test data
  } // topology
} // pylith

// RefineUniform ---------------------------------------------------------------
class pylith::topology::TestRefineUniform : public CppUnit::TestFixture
{ // class TestRefineUniform

  // CPPUNIT TEST SUITE /////////////////////////////////////////////////
  CPPUNIT_TEST_SUITE( TestRefineUniform );

  CPPUNIT_TEST( testConstructor );

  CPPUNIT_TEST( testRefineTri3Level2 );
  CPPUNIT_TEST( testRefineTri3Level2Fault1 );

  CPPUNIT_TEST( testRefineQuad4Level2 );
  CPPUNIT_TEST( testRefineQuad4Level2Fault1 );

  CPPUNIT_TEST( testRefineTet4Level2 );
  CPPUNIT_TEST( testRefineTet4Level2Fault1 );

  CPPUNIT_TEST( testRefineHex8Level2 );
  CPPUNIT_TEST( testRefineHex8Level2Fault1 );

  CPPUNIT_TEST_SUITE_END();

  // PUBLIC METHODS /////////////////////////////////////////////////////
public :

  /// Test constructor.
  void testConstructor(void);

  /// Test refine() with level 2, tri3 cells, and no fault.
  void testRefineTri3Level2(void);

  /// Test refine() with level 2, tri3 cells, and one fault.
  void testRefineTri3Level2Fault1(void);

  /// Test refine() with level 2, quad4 cells, and no fault.
  void testRefineQuad4Level2(void);

  /// Test refine() with level 2, quad4 cells, and one fault.
  void testRefineQuad4Level2Fault1(void);

  /// Test refine() with level 2, tet4 cells, and no fault.
  void testRefineTet4Level2(void);

  /// Test refine() with level 2, tet4 cells, and one fault.
  void testRefineTet4Level2Fault1(void);

  /// Test refine() with level 2, hex8 cells, and no fault.
  void testRefineHex8Level2(void);

  /// Test refine() with level 2, hex8 cells, and one fault.
  void testRefineHex8Level2Fault1(void);

// PRIVATE METHODS //////////////////////////////////////////////////////
private :

  /** Setup mesh.
   *
   * @mesh Mesh to setup.
   * @param data Test data.
   */
  void _setupMesh(Mesh* const mesh,
		  const MeshDataCohesive& data);

  /** Test refine().
   *
   * @param data Test data.
   */
  void _testRefine(const MeshDataCohesive& data);

}; // class TestRefineUniform

#endif // pylith_topology_testrefineuniform_hh


// End of file 
