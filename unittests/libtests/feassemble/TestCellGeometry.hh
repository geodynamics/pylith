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
 * @file unittests/libtests/feassemble/TestCellGeometry.hh
 *
 * @brief C++ TestCellGeometry object
 *
 * C++ unit testing for CellGeometry.
 */

#if !defined(pylith_feassemble_testcellgeometry_hh)
#define pylith_feassemble_testcellgeometry_hh

#include <cppunit/extensions/HelperMacros.h>

/// Namespace for pylith package
namespace pylith {
  namespace feassemble {
    class TestCellGeometry;

    class CellGeometry;
    class CellGeomData;
  } // feassemble
} // pylith

/// C++ unit testing for TestCellGeometry
class pylith::feassemble::TestCellGeometry : public CppUnit::TestFixture
{ // class TestCellGeometry

  // CPPUNIT TEST SUITE /////////////////////////////////////////////////
  CPPUNIT_TEST_SUITE( TestCellGeometry );

  CPPUNIT_TEST( testOrient1D );
  CPPUNIT_TEST( testOrient2D );
  CPPUNIT_TEST( testOrient3D );

  CPPUNIT_TEST_SUITE_END();

// PUBLIC METHODS ///////////////////////////////////////////////////////
public :

  /// Setup data.
  void setUp(void);

  /// Tear down data.
  void tearDown(void);

  /// Test _orient1D().
  void testOrient1D(void);

  /// Test _orient2D().
  void testOrient2D(void);

  /// Test _orient3D().
  void testOrient3D(void);

  /// Test clone().
  void testClone(void);

  /// Test cellDim().
  void testCellDim(void);

  /// Test spaceDim().
  void testSpaceDim(void);

  /// Test numCorners().
  void testNumCorners(void);

  /// Test orientFn.
  void testOrientFn(void);

  /// Test jacobian().
  void testJacobian(void);

// PROTECTED MEMBERS ////////////////////////////////////////////////////
protected :

  CellGeometry* _object; ///< Test subject.
  CellGeomData* _data; ///< Test data.

}; // class TestCellGeometry

#endif // pylith_feassemble_testcellgeometry_hh

// End of file 
