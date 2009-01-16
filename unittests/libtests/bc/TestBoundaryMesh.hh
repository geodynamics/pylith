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
 * @file unittests/libtests/bc/TestBoundaryMesh.hh
 *
 * @brief C++ TestBoundaryMesh object.
 *
 * C++ unit testing for BoundaryMesh.
 */

#if !defined(pylith_bc_testboundarymesh_hh)
#define pylith_bc_testboundarymesh_hh

#include <cppunit/extensions/HelperMacros.h>

/// Namespace for pylith package
namespace pylith {
  namespace bc {
    class TestBoundaryMesh;

    class BoundaryMeshData;
  } // bc
} // pylith

/// C++ unit testing for BoundaryMesh.
class pylith::bc::TestBoundaryMesh : public CppUnit::TestFixture
{ // class TestBoundaryMesh

  // PUBLIC METHODS /////////////////////////////////////////////////////
public :

  /// Setup testing data.
  void setUp(void);

  /// Tear down testing data.
  void tearDown(void);

  /// Test submesh() without fault.
  void testSubmesh(void);

  /// Test submesh() with fault().
  void testSubmeshFault(void);

  // PROTECTED MEMBERS //////////////////////////////////////////////////
protected :

  BoundaryMeshData* _data; ///< Data for testing
  bool _flipFault; ///< If true, flip fault orientation.

}; // class TestBoundaryMesh

#endif // pylith_bc_boundarymesh_hh


// End of file 
