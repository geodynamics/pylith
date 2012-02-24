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
// Copyright (c) 2010-2012 University of California, Davis
//
// See COPYING for license information.
//
// ----------------------------------------------------------------------
//

/**
 * @file unittests/libtests/bc/TestNeumannHex8.hh
 *
 * @brief C++ TestNeumann object.
 *
 * C++ unit testing for Neumann for mesh with 3-D hex cells.
 */

#if !defined(pylith_bc_testneumannhex8_hh)
#define pylith_bc_testneumannhex8_hh

#include "TestNeumann.hh" // ISA TestNeumann

/// Namespace for pylith package
namespace pylith {
  namespace bc {
    class TestNeumannHex8;
  } // bc
} // pylith

/// C++ unit testing for Neumann for mesh with 3-D hex cells.
class pylith::bc::TestNeumannHex8 : public TestNeumann
{ // class TestNeumann

  // CPPUNIT TEST SUITE /////////////////////////////////////////////////
  CPPUNIT_TEST_SUITE( TestNeumannHex8 );

  CPPUNIT_TEST( testInitialize );
  CPPUNIT_TEST( testIntegrateResidual );

  CPPUNIT_TEST_SUITE_END();

  // PUBLIC METHODS /////////////////////////////////////////////////////
public :

  /// Setup testing data.
  void setUp(void);

}; // class TestNeumannHex8

#endif // pylith_bc_neumannhex8_hh


// End of file 
