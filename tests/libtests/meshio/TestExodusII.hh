// -*- C++ -*-
//
// ----------------------------------------------------------------------
//
// Brad T. Aagaard, U.S. Geological Survey
// Charles A. Williams, GNS Science
// Matthew G. Knepley, University at Buffalo
//
// This code was developed as part of the Computational Infrastructure
// for Geodynamics (http://geodynamics.org).
//
// Copyright (c) 2010-2021 University of California, Davis
//
// See LICENSE.md for license information.
//
// ----------------------------------------------------------------------
//

/**
 * @file tests/libtests/meshio/TestExodusII.hh
 *
 * @brief C++ TestExodusII object
 *
 * C++ unit testing for ExodusII.
 */

#if !defined(pylith_meshio_testexodusii_hh)
#define pylith_meshio_testexodusii_hh

// Include directives ---------------------------------------------------
#include <cppunit/extensions/HelperMacros.h>

// Forward declarations -------------------------------------------------
namespace pylith {
  namespace meshio {
    class TestExodusII;
  } // meshio
} // pylith

// TestExodusII ------------------------------------------------------
class pylith::meshio::TestExodusII : public CppUnit::TestFixture
{ // class TestExodusII

  // CPPUNIT TEST SUITE /////////////////////////////////////////////////
  CPPUNIT_TEST_SUITE( TestExodusII );

  CPPUNIT_TEST( testConstructor );
  CPPUNIT_TEST( testFilename );
  CPPUNIT_TEST( testOpenClose );
  CPPUNIT_TEST( testHasDim );
  CPPUNIT_TEST( testHasAtt );
  CPPUNIT_TEST( testHasVar );
  CPPUNIT_TEST( testGetVarDouble );
  CPPUNIT_TEST( testGetVarInt );
  CPPUNIT_TEST( testGetVarString );

  CPPUNIT_TEST_SUITE_END();

  // PUBLIC METHODS /////////////////////////////////////////////////////
public :

  /// Test constructor
  void testConstructor(void);

  /// Test filename()
  void testFilename(void);

  /// Test open() and close().
  void testOpenClose(void);

  /// Test hasDim()
  void testHasDim(void);

  /// Test hasAtt()
  void testHasAtt(void);

  /// Test hasVar()
  void testHasVar(void);

  /// Test getVar(PylithScalar*)
  void testGetVarDouble(void);

  /// Test getVar(int*)
  void testGetVarInt(void);

  /// Test getVar(string_vector)
  void testGetVarString(void);

}; // class TestExodusII

#endif // pylith_meshio_testexodusii_hh


// End of file 
