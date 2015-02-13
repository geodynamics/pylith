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
 * @file unittests/libtests/bc/TestBoundaryMeshTri3.hh
 *
 * @brief C++ TestBoundaryMesh object.
 *
 * Test cases for C++ unit testing of submesh().
 */

#if !defined(pylith_bc_testboundarymeshcases_hh)
#define pylith_bc_testboundarymeshcases_hh

#include "TestBoundaryMesh.hh" // ISA TestBoundaryMesh

/// Namespace for pylith package
namespace pylith {
  namespace bc {
    class TestBoundaryMeshTri3;
    class TestBoundaryMeshQuad4;
    class TestBoundaryMeshTet4;
    class TestBoundaryMeshHex8;
  } // bc
} // pylith

// ----------------------------------------------------------------------
/// C++ unit testing of submesh() for mesh with 2-D tri cells.
class pylith::bc::TestBoundaryMeshTri3 : public TestBoundaryMesh
{ // class TestBoundaryMesh

  // CPPUNIT TEST SUITE /////////////////////////////////////////////////
  CPPUNIT_TEST_SUITE(TestBoundaryMeshTri3);

  CPPUNIT_TEST( testSubmesh );
  CPPUNIT_TEST( testSubmeshFault );

  CPPUNIT_TEST_SUITE_END();

  // PUBLIC METHODS /////////////////////////////////////////////////////
public :

  /// Setup testing data.
  void setUp(void);

}; // class TestBoundaryMeshTri3


// ----------------------------------------------------------------------
/// C++ unit testing of submesh() for mesh with 2-D quad cells.
class pylith::bc::TestBoundaryMeshQuad4 : public TestBoundaryMesh
{ // class TestBoundaryMesh

  // CPPUNIT TEST SUITE /////////////////////////////////////////////////
  CPPUNIT_TEST_SUITE(TestBoundaryMeshQuad4);

  CPPUNIT_TEST( testSubmesh );
  CPPUNIT_TEST( testSubmeshFault );

  CPPUNIT_TEST_SUITE_END();

  // PUBLIC METHODS /////////////////////////////////////////////////////
public :

  /// Setup testing data.
  void setUp(void);

}; // class TestBoundaryMeshQuad4


// ----------------------------------------------------------------------
/// C++ unit testing of submesh() for mesh with 3-D tet cells.
class pylith::bc::TestBoundaryMeshTet4 : public TestBoundaryMesh
{ // class TestBoundaryMesh

  // CPPUNIT TEST SUITE /////////////////////////////////////////////////
  CPPUNIT_TEST_SUITE(TestBoundaryMeshTet4);

  CPPUNIT_TEST( testSubmesh );
  CPPUNIT_TEST( testSubmeshFault );

  CPPUNIT_TEST_SUITE_END();

  // PUBLIC METHODS /////////////////////////////////////////////////////
public :

  /// Setup testing data.
  void setUp(void);

}; // class TestBoundaryMeshTet4


// ----------------------------------------------------------------------
/// C++ unit testing of submesh() for mesh with 3-D hex cells.
class pylith::bc::TestBoundaryMeshHex8 : public TestBoundaryMesh
{ // class TestBoundaryMesh

  // CPPUNIT TEST SUITE /////////////////////////////////////////////////
  CPPUNIT_TEST_SUITE(TestBoundaryMeshHex8);

  CPPUNIT_TEST( testSubmesh );
  CPPUNIT_TEST( testSubmeshFault );

  CPPUNIT_TEST_SUITE_END();

  // PUBLIC METHODS /////////////////////////////////////////////////////
public :

  /// Setup testing data.
  void setUp(void);

}; // class TestBoundaryMeshHex8


#endif // pylith_bc_boundarymeshcases_hh


// End of file 
