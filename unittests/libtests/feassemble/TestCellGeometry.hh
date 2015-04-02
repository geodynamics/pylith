// -*- C++ -*-
//
// ----------------------------------------------------------------------
//
// Brad T. Aagaard, U.S. Geological Survey
// Charles A. Williams, GNS Science
// Matthew G. Knepley, University of Chicago
//
// This code was developed as part of the Computational Infrastructure
// for Geodynamics (http://geodynamics.org).
//
// Copyright (c) 2010-2015 University of California, Davis
//
// See COPYING for license information.
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

#include "pylith/feassemble/feassemblefwd.hh" // forward declarations

/// Namespace for pylith package
namespace pylith {
  namespace feassemble {
    class TestCellGeometry;

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
