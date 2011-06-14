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
// Copyright (c) 2010-2011 University of California, Davis
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
 * C++ unit testing of submesh() for mesh with 2-D tri cells.
 */

#if !defined(pylith_bc_testboundarymeshtri3_hh)
#define pylith_bc_testboundarymeshtri3_hh

#include "TestBoundaryMesh.hh" // ISA TestBoundaryMesh

/// Namespace for pylith package
namespace pylith {
  namespace bc {
    class TestBoundaryMeshTri3;
  } // bc
} // pylith

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

#endif // pylith_bc_boundarymeshtri3_hh


// End of file 
