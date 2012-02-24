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
 * @file unittests/libtests/bc/TestBoundaryMeshTet4.hh
 *
 * @brief C++ TestBoundaryMesh object.
 *
 * C++ unit testing of submesh() for mesh with 2-D tri cells.
 */

#if !defined(pylith_bc_testboundarymeshtet4_hh)
#define pylith_bc_testboundarymeshtet4_hh

#include "TestBoundaryMesh.hh" // ISA TestBoundaryMesh

/// Namespace for pylith package
namespace pylith {
  namespace bc {
    class TestBoundaryMeshTet4;
  } // bc
} // pylith

/// C++ unit testing of submesh() for mesh with 2-D tri cells.
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

#endif // pylith_bc_boundarymeshtet4_hh


// End of file 
