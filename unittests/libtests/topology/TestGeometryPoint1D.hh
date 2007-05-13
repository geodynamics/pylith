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
 * @file unittests/libtests/topology/TestGeometryPoint1D.hh
 *
 * @brief C++ TestGeometryPoint1D object
 *
 * C++ unit testing for TopologyAscii.
 */

#if !defined(pylith_topology_testgeometrypoint1d_hh)
#define pylith_topology_testgeometrypoint1d_hh

#include "TestCellGeometry.hh"

/// Namespace for pylith package
namespace pylith {
  namespace topology {
    class TestGeometryPoint1D;
    class CellGeomData;
  } // topology
} // pylith

/// C++ unit testing for Quadrature1D
class pylith::topology::TestGeometryPoint1D : public TestCellGeometry
{ // class TestGeometryPoint1D

  // CPPUNIT TEST SUITE /////////////////////////////////////////////////
  CPPUNIT_TEST_SUITE( TestGeometryPoint1D );
  CPPUNIT_TEST( testConstructor );
  CPPUNIT_TEST( testCellDim );
  CPPUNIT_TEST( testSpaceDim );
  CPPUNIT_TEST( testNumCorners );
  CPPUNIT_TEST( testGeomLowerDim );
  CPPUNIT_TEST( testJacobian );
  CPPUNIT_TEST_SUITE_END();

  // PUBLIC METHODS /////////////////////////////////////////////////////
public :

  /// Test constructor
  void testConstructor(void);

  /// Test cellDim()
  void testCellDim(void);

  /// Test spaceDim()
  void testSpaceDim(void);

  /// Test numCorners()
  void testNumCorners(void);

  /// Test geometryLowerDim().
  void testGeomLowerDim(void);

  /// Test jacobian().
  void testJacobian(void);

}; // class TestGeometryPoint1D

#endif // pylith_topology_testgeometrypoint1d_hh

// End of file 
