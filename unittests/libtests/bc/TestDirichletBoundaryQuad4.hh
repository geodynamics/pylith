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
 * @file unittests/libtests/bc/TestDirichletBoundaryQuad4.hh
 *
 * @brief C++ TestDirichletBoundary object.
 *
 * C++ unit testing for DirichletBoundary for mesh with 1-D line cells.
 */

#if !defined(pylith_bc_testdirichletboundaryquad4_hh)
#define pylith_bc_testdirichletboundaryquad4_hh

#include "TestDirichletBoundary.hh" // ISA TestDirichletBoundary

/// Namespace for pylith package
namespace pylith {
  namespace bc {
    class TestDirichletBoundaryQuad4;
  } // bc
} // pylith

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

#endif // pylith_bc_dirichletboundaryquad4_hh


// End of file 
