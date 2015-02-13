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
 * @file unittests/libtests/bc/TestNeumannCases.hh
 *
 * @brief C++ TestNeumann object.
 *
 * Test cases for C++ unit testing for Neumann.
 */

#if !defined(pylith_bc_testneumanncases_hh)
#define pylith_bc_testneumanncases_hh

#include "TestNeumann.hh" // ISA TestNeumann

/// Namespace for pylith package
namespace pylith {
  namespace bc {
    class TestNeumannTri3;
    class TestNeumannQuad4;
    class TestNeumannTet4;
    class TestNeumannHex8;
  } // bc
} // pylith

// ----------------------------------------------------------------------
/// C++ unit testing for Neumann for mesh with 2-D tri cells.
class pylith::bc::TestNeumannTri3 : public TestNeumann
{ // class TestNeumann

  // CPPUNIT TEST SUITE /////////////////////////////////////////////////
  CPPUNIT_TEST_SUITE( TestNeumannTri3 );

  CPPUNIT_TEST( testInitialize );
  CPPUNIT_TEST( testIntegrateResidual );

  CPPUNIT_TEST_SUITE_END();

  // PUBLIC METHODS /////////////////////////////////////////////////////
public :

  /// Setup testing data.
  void setUp(void);

}; // class TestNeumannTri3


// ----------------------------------------------------------------------
/// C++ unit testing for Neumann for mesh with 2-D quad cells.
class pylith::bc::TestNeumannQuad4 : public TestNeumann
{ // class TestNeumann

  // CPPUNIT TEST SUITE /////////////////////////////////////////////////
  CPPUNIT_TEST_SUITE( TestNeumannQuad4 );

  CPPUNIT_TEST( testInitialize );
  CPPUNIT_TEST( testIntegrateResidual );

  CPPUNIT_TEST_SUITE_END();

  // PUBLIC METHODS /////////////////////////////////////////////////////
public :

  /// Setup testing data.
  void setUp(void);

}; // class TestNeumannQuad4


// ----------------------------------------------------------------------
/// C++ unit testing for Neumann for mesh with 3-D tet cells.
class pylith::bc::TestNeumannTet4 : public TestNeumann
{ // class TestNeumann

  // CPPUNIT TEST SUITE /////////////////////////////////////////////////
  CPPUNIT_TEST_SUITE( TestNeumannTet4 );

  CPPUNIT_TEST( testInitialize );
  CPPUNIT_TEST( testIntegrateResidual );

  CPPUNIT_TEST_SUITE_END();

  // PUBLIC METHODS /////////////////////////////////////////////////////
public :

  /// Setup testing data.
  void setUp(void);

}; // class TestNeumannTet4


// ----------------------------------------------------------------------
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


#endif // pylith_bc_neumanncases_hh


// End of file 
