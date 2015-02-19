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

  CPPUNIT_TEST( testRefineTri3Level1 );
  CPPUNIT_TEST( testRefineTri3Level1Fault1 );

  CPPUNIT_TEST( testRefineQuad4Level1 );
  CPPUNIT_TEST( testRefineQuad4Level1Fault1 );

  CPPUNIT_TEST( testRefineTet4Level1 );
  CPPUNIT_TEST( testRefineTet4Level1Fault1 );

  CPPUNIT_TEST( testRefineHex8Level1 );
  CPPUNIT_TEST( testRefineHex8Level1Fault1 );

  CPPUNIT_TEST_SUITE_END();

  // PUBLIC METHODS /////////////////////////////////////////////////////
public :

  /// Test constructor.
  void testConstructor(void);

  /// Test refine() with level 1, tri3 cells, and no fault.
  void testRefineTri3Level1(void);

  /// Test refine() with level 1, tri3 cells, and one fault.
  void testRefineTri3Level1Fault1(void);

  /// Test refine() with level 1, quad4 cells, and no fault.
  void testRefineQuad4Level1(void);

  /// Test refine() with level 1, quad4 cells, and one fault.
  void testRefineQuad4Level1Fault1(void);

  /// Test refine() with level 1, tet4 cells, and no fault.
  void testRefineTet4Level1(void);

  /// Test refine() with level 1, tet4 cells, and one fault.
  void testRefineTet4Level1Fault1(void);

  /// Test refine() with level 1, hex8 cells, and no fault.
  void testRefineHex8Level1(void);

  /// Test refine() with level 1, hex8 cells, and one fault.
  void testRefineHex8Level1Fault1(void);

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
   * @param isSimplexMesh True if mesh is simplex (tri,tet), false if not (quad,hex).
   */
  void _testRefine(const MeshDataCohesive& data,
		   const bool isSimplexMesh);

}; // class TestRefineUniform

#endif // pylith_topology_testrefineuniform_hh


// End of file 
