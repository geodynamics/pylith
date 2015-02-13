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
 * @file unittests/libtests/topology/TestReverseCuthillMcKee.hh
 *
 * @brief C++ TestReverseCuthillMcKee object
 *
 * C++ unit testing for ReverseCuthillMcKee.
 */

#if !defined(pylith_topology_testreversecuthillmckee_hh)
#define pylith_topology_testreversecuthillmckee_hh

// Include directives ---------------------------------------------------
#include <cppunit/extensions/HelperMacros.h>

#include "pylith/topology/topologyfwd.hh" // USES Mesh

// Forward declarations -------------------------------------------------
/// Namespace for pylith package
namespace pylith {
  namespace topology {
    class TestReverseCuthillMcKee;
  } // topology
} // pylith

// ReverseCuthillMcKee ---------------------------------------------------------------
class pylith::topology::TestReverseCuthillMcKee : public CppUnit::TestFixture
{ // class TestReverseCuthillMcKee

  // CPPUNIT TEST SUITE /////////////////////////////////////////////////
  CPPUNIT_TEST_SUITE( TestReverseCuthillMcKee );

  CPPUNIT_TEST( testReorderTri3 );
  CPPUNIT_TEST( testReorderTri3Fault );

  CPPUNIT_TEST( testReorderQuad4 );
  CPPUNIT_TEST( testReorderQuad4Fault );

  CPPUNIT_TEST( testReorderTet4 );
  CPPUNIT_TEST( testReorderTet4Fault );

  CPPUNIT_TEST( testReorderHex8 );
  CPPUNIT_TEST( testReorderHex8Fault );

  CPPUNIT_TEST_SUITE_END();

  // PUBLIC METHODS /////////////////////////////////////////////////////
public :

  /// Test reorder() with tri3 cells and no fault.
  void testReorderTri3(void);

  /// Test reorder() with tri3 cells and one fault.
  void testReorderTri3Fault(void);

  /// Test reorder() with quad4 cells and no fault.
  void testReorderQuad4(void);

  /// Test reorder() with quad4 cells and one fault.
  void testReorderQuad4Fault(void);

  /// Test reorder() with tet4 cells and no fault.
  void testReorderTet4(void);

  /// Test reorder() with tet4 cells and one fault.
  void testReorderTet4Fault(void);

  /// Test reorder() with hex8 cells and no fault.
  void testReorderHex8(void);

  /// Test reorder() with hex8 cells and one fault.
  void testReorderHex8Fault(void);

// PRIVATE METHODS //////////////////////////////////////////////////////
private :

  /** Setup mesh.
   *
   * @mesh Mesh to setup.
   * @param filename Mesh filename.
   * @param faultGroup Name of fault group.
   */
  void _setupMesh(Mesh* const mesh,
		  const char* filename,
		  const char* faultGroup =0);

  /** Test reorder().
   *
   * @param filename Mesh filename.
   * @param faultGroup Name of fault group.
   */
  void _testReorder(const char* filename,
		    const char* faultGroup =0);

}; // class TestReverseCuthillMcKee

#endif // pylith_topology_testreversecuthillmckee_hh


// End of file 
