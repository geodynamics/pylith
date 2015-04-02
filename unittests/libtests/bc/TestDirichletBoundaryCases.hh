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
 * @file unittests/libtests/bc/TestDirichletBoundaryCases.hh
 *
 * @brief C++ TestDirichletBoundary object.
 *
 * Test cases for C++ unit testing for DirichletBoundary for mesh.
 */

#if !defined(pylith_bc_testdirichletboundarycases_hh)
#define pylith_bc_testdirichletboundarycases_hh

#include "TestDirichletBoundary.hh" // ISA TestDirichletBoundary

/// Namespace for pylith package
namespace pylith {
  namespace bc {
    class TestDirichletBoundaryTri3;
    class TestDirichletBoundaryQuad4;
    class TestDirichletBoundaryTet4;
    class TestDirichletBoundaryHex8;
  } // bc
} // pylith

// ----------------------------------------------------------------------
/// C++ unit testing for DirichletBoundary for mesh with 2-D tri cells.
class pylith::bc::TestDirichletBoundaryTri3 : public TestDirichletBoundary
{ // class TestDirichletBoundary

  // CPPUNIT TEST SUITE /////////////////////////////////////////////////
  CPPUNIT_TEST_SUB_SUITE( TestDirichletBoundaryTri3, TestDirichletBoundary );
  CPPUNIT_TEST( testInitialize );
  CPPUNIT_TEST_SUITE_END();

  // PUBLIC METHODS /////////////////////////////////////////////////////
public :

  /// Setup testing data.
  void setUp(void);

}; // class TestDirichletBoundaryTri3


// ----------------------------------------------------------------------
/// C++ unit testing for DirichletBoundary for mesh with 2-D quad cells.
class pylith::bc::TestDirichletBoundaryQuad4 : public TestDirichletBoundary
{ // class TestDirichletBoundary

  // CPPUNIT TEST SUITE /////////////////////////////////////////////////
  CPPUNIT_TEST_SUB_SUITE( TestDirichletBoundaryQuad4, TestDirichletBoundary );
  CPPUNIT_TEST( testInitialize );
  CPPUNIT_TEST_SUITE_END();

  // PUBLIC METHODS /////////////////////////////////////////////////////
public :

  /// Setup testing data.
  void setUp(void);

}; // class TestDirichletBoundaryQuad4


// ----------------------------------------------------------------------
/// C++ unit testing for DirichletBoundary for mesh with 3-D tet cells.
class pylith::bc::TestDirichletBoundaryTet4 : public TestDirichletBoundary
{ // class TestDirichletBoundary

  // CPPUNIT TEST SUITE /////////////////////////////////////////////////
  CPPUNIT_TEST_SUB_SUITE( TestDirichletBoundaryTet4, TestDirichletBoundary );
  CPPUNIT_TEST( testInitialize );
  CPPUNIT_TEST_SUITE_END();

  // PUBLIC METHODS /////////////////////////////////////////////////////
public :

  /// Setup testing data.
  void setUp(void);

}; // class TestDirichletBoundaryTet4


// ----------------------------------------------------------------------
/// C++ unit testing for DirichletBoundary for mesh with 3-D hex cells.
class pylith::bc::TestDirichletBoundaryHex8 : public TestDirichletBoundary
{ // class TestDirichletBoundary

  // CPPUNIT TEST SUITE /////////////////////////////////////////////////
  CPPUNIT_TEST_SUB_SUITE( TestDirichletBoundaryHex8, TestDirichletBoundary );
  CPPUNIT_TEST( testInitialize );
  CPPUNIT_TEST_SUITE_END();

  // PUBLIC METHODS /////////////////////////////////////////////////////
public :

  /// Setup testing data.
  void setUp(void);

}; // class TestDirichletBoundaryHex8


#endif // pylith_bc_dirichletboundarycases_hh


// End of file 
