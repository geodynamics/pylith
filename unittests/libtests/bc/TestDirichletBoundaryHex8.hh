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
// Copyright (c) 2010-2011 University of California, Davis
//
// See COPYING for license information.
//
// ----------------------------------------------------------------------
//

/**
 * @file unittests/libtests/bc/TestDirichletBoundaryHex8.hh
 *
 * @brief C++ TestDirichletBoundary object.
 *
 * C++ unit testing for DirichletBoundary for mesh with 1-D line cells.
 */

#if !defined(pylith_bc_testdirichletboundaryhex8_hh)
#define pylith_bc_testdirichletboundaryhex8_hh

#include "TestDirichletBoundary.hh" // ISA TestDirichletBoundary

/// Namespace for pylith package
namespace pylith {
  namespace bc {
    class TestDirichletBoundaryHex8;
  } // bc
} // pylith

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

#endif // pylith_bc_dirichletboundaryhex8_hh


// End of file 
