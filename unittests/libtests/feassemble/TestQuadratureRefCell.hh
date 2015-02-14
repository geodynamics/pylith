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
 * @file unittests/libtests/feassemble/TestQuadratureRefCell.hh
 *
 * @brief C++ TestQuadratureRefCell object
 *
 * C++ unit testing for QuadratureRefCell.
 */

#if !defined(pylith_feassemble_testquadraturerefcell_hh)
#define pylith_feassemble_testquadraturerefcell_hh

#include <cppunit/extensions/HelperMacros.h>

/// Namespace for pylith package
namespace pylith {
  namespace feassemble {
    class TestQuadratureRefCell;
    class QuadratureRefCellData;
  } // feassemble
} // pylith

/// C++ unit testing for QuadratureRefCell
class pylith::feassemble::TestQuadratureRefCell : public CppUnit::TestFixture
{ // class TestQuadratureRefCell

  // CPPUNIT TEST SUITE /////////////////////////////////////////////////
  CPPUNIT_TEST_SUITE( TestQuadratureRefCell );

  CPPUNIT_TEST( testConstructor );
  CPPUNIT_TEST( testMinJacobian );
  CPPUNIT_TEST( testRefGeometry );
  CPPUNIT_TEST( testInitialize );

  CPPUNIT_TEST_SUITE_END();

  // PUBLIC METHODS /////////////////////////////////////////////////////
public :

  /// Test constructor.
  void testConstructor(void);

  /// Test minJacobian()
  void testMinJacobian(void);

  /// Test refGeometry()
  void testRefGeometry(void);

  /// Test initialize()
  void testInitialize(void);

}; // class TestQuadratureRefCell

#endif // pylith_feassemble_testquadraturerefcell_hh

// End of file 
