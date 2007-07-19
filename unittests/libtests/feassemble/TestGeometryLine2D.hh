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
 * @file unittests/libtests/feassemble/TestGeometryLine2D.hh
 *
 * @brief C++ TestGeometryLine2D object
 *
 * C++ unit testing for TopologyAscii.
 */

#if !defined(pylith_feassemble_testgeometryline2d_hh)
#define pylith_feassemble_testgeometryline2d_hh

#include "TestCellGeometry.hh"

/// Namespace for pylith package
namespace pylith {
  namespace feassemble {
    class TestGeometryLine2D;
  } // feassemble
} // pylith

/// C++ unit testing for GeometryLine2D
class pylith::feassemble::TestGeometryLine2D : public TestCellGeometry
{ // class TestGeometryLine2D

  // CPPUNIT TEST SUITE /////////////////////////////////////////////////
  CPPUNIT_TEST_SUITE( TestGeometryLine2D );

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

}; // class TestGeometryLine2D

#endif // pylith_feassemble_testgeometryline2d_hh

// End of file 
