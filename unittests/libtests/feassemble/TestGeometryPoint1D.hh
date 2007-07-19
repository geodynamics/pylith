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
 * @file unittests/libtests/feassemble/TestGeometryPoint1D.hh
 *
 * @brief C++ TestGeometryPoint1D object
 *
 * C++ unit testing for TopologyAscii.
 */

#if !defined(pylith_feassemble_testgeometrypoint1d_hh)
#define pylith_feassemble_testgeometrypoint1d_hh

#include "TestCellGeometry.hh"

/// Namespace for pylith package
namespace pylith {
  namespace feassemble {
    class TestGeometryPoint1D;
  } // feassemble
} // pylith

/// C++ unit testing for GeometryPoint1D
class pylith::feassemble::TestGeometryPoint1D : public TestCellGeometry
{ // class TestGeometryPoint1D

  // CPPUNIT TEST SUITE /////////////////////////////////////////////////
  CPPUNIT_TEST_SUITE( TestGeometryPoint1D );

  CPPUNIT_TEST( testConstructor );
  CPPUNIT_TEST( testClone );
  CPPUNIT_TEST( testGeomLowerDim );

  CPPUNIT_TEST( testCellDim );
  CPPUNIT_TEST( testSpaceDim );
  CPPUNIT_TEST( testNumCorners );
  CPPUNIT_TEST( testOrientFn );
  CPPUNIT_TEST( testJacobian );

  CPPUNIT_TEST_SUITE_END();

  // PUBLIC METHODS /////////////////////////////////////////////////////
public :

  /// Setup data.
  void setUp(void);

  /// Test constructor.
  void testConstructor(void);

  /// Test geometryLowerDim().
  void testGeomLowerDim(void);

  /// Test jacobian().
  void testJacobian(void);

}; // class TestGeometryPoint1D

#endif // pylith_feassemble_testgeometrypoint1d_hh

// End of file 
