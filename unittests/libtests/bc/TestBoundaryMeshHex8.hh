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
 * @file unittests/libtests/bc/TestBoundaryMeshHex8.hh
 *
 * @brief C++ TestBoundaryMesh object.
 *
 * C++ unit testing of submesh() for mesh with 2-D tri cells.
 */

#if !defined(pylith_bc_testboundarymeshhex8_hh)
#define pylith_bc_testboundarymeshhex8_hh

#include "TestBoundaryMesh.hh" // ISA TestBoundaryMesh

/// Namespace for pylith package
namespace pylith {
  namespace bc {
    class TestBoundaryMeshHex8;
  } // bc
} // pylith

/// C++ unit testing of submesh() for mesh with 2-D tri cells.
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

#endif // pylith_bc_boundarymeshhex8_hh


// End of file 
