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
 * @file unittests/libtests/feassemble/TestGeometryHex3D.hh
 *
 * @brief C++ TestGeometryHex3D object
 *
 * C++ unit testing for TopologyAscii.
 */

#if !defined(pylith_feassemble_testgeometryhex3d_hh)
#define pylith_feassemble_testgeometryhex3d_hh

#include "TestCellGeometry.hh"

/// Namespace for pylith package
namespace pylith {
  namespace feassemble {
    class TestGeometryHex3D;
  } // feassemble
} // pylith

/// C++ unit testing for GeometryHex3D
class pylith::feassemble::TestGeometryHex3D : public TestCellGeometry
{ // class TestGeometryHex3D

  // CPPUNIT TEST SUITE /////////////////////////////////////////////////
  CPPUNIT_TEST_SUITE( TestGeometryHex3D );

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

}; // class TestGeometryHex3D

#endif // pylith_feassemble_testgeometryhex3d_hh

// End of file 
