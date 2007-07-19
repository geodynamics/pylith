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
 * @file unittests/libtests/feassemble/TestGeometryLine3D.hh
 *
 * @brief C++ TestGeometryLine3D object
 *
 * C++ unit testing for TopologyAscii.
 */

#if !defined(pylith_feassemble_testgeometryline3d_hh)
#define pylith_feassemble_testgeometryline3d_hh

#include "TestCellGeometry.hh"

/// Namespace for pylith package
namespace pylith {
  namespace feassemble {
    class TestGeometryLine3D;
  } // feassemble
} // pylith

/// C++ unit testing for GeometryLine3D
class pylith::feassemble::TestGeometryLine3D : public TestCellGeometry
{ // class TestGeometryLine3D

  // CPPUNIT TEST SUITE /////////////////////////////////////////////////
  CPPUNIT_TEST_SUITE( TestGeometryLine3D );

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

}; // class TestGeometryLine3D

#endif // pylith_feassemble_testgeometryline3d_hh

// End of file 
