// -*- C++ -*-
//
// ----------------------------------------------------------------------
//
//                           Brad T. Aagaard
//                        U.S. Geological Survey
//
// {LicenseText}
//
// ----------------------------------------------------------------------
//

/**
 * @file unittests/libtests/bc/TestBoundaryMeshQuad4.hh
 *
 * @brief C++ TestBoundaryMesh object.
 *
 * C++ unit testing of submesh() for mesh with 2-D tri cells.
 */

#if !defined(pylith_bc_testboundarymeshquad4_hh)
#define pylith_bc_testboundarymeshquad4_hh

#include "TestBoundaryMesh.hh" // ISA TestBoundaryMesh

/// Namespace for pylith package
namespace pylith {
  namespace bc {
    class TestBoundaryMeshQuad4;
  } // bc
} // pylith

/// C++ unit testing of submesh() for mesh with 2-D tri cells.
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

#endif // pylith_bc_boundarymeshquad4_hh


// End of file 
