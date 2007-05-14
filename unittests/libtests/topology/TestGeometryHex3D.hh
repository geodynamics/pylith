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
 * @file unittests/libtests/topology/TestGeometryHex3D.hh
 *
 * @brief C++ TestGeometryHex3D object
 *
 * C++ unit testing for TopologyAscii.
 */

#if !defined(pylith_topology_testgeometryhex3d_hh)
#define pylith_topology_testgeometryhex3d_hh

#include "TestCellGeometry.hh"

/// Namespace for pylith package
namespace pylith {
  namespace topology {
    class TestGeometryHex3D;
    class CellGeomData;
  } // topology
} // pylith

/// C++ unit testing for Quadrature3D
class pylith::topology::TestGeometryHex3D : public TestCellGeometry
{ // class TestGeometryHex3D

  // CPPUNIT TEST SUITE /////////////////////////////////////////////////
  CPPUNIT_TEST_SUITE( TestGeometryHex3D );
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

}; // class TestGeometryHex3D

#endif // pylith_topology_testgeometryhex3d_hh

// End of file 
