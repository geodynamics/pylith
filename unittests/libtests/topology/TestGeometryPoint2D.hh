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
 * @file unittests/libtests/topology/TestGeometryPoint2D.hh
 *
 * @brief C++ TestGeometryPoint2D object
 *
 * C++ unit testing for TopologyAscii.
 */

#if !defined(pylith_topology_testgeometrypoint2d_hh)
#define pylith_topology_testgeometrypoint2d_hh

#include "TestCellGeometry.hh"

/// Namespace for pylith package
namespace pylith {
  namespace topology {
    class TestGeometryPoint2D;
    class CellGeomData;
  } // topology
} // pylith

/// C++ unit testing for Quadrature2D
class pylith::topology::TestGeometryPoint2D : public TestCellGeometry
{ // class TestGeometryPoint2D

  // CPPUNIT TEST SUITE /////////////////////////////////////////////////
  CPPUNIT_TEST_SUITE( TestGeometryPoint2D );
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

}; // class TestGeometryPoint2D

#endif // pylith_topology_testgeometrypoint2d_hh

// End of file 
