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
 * @file unittests/libtests/bc/TestDirichletBCTet4.hh
 *
 * @brief C++ TestDirichletBC object.
 *
 * C++ unit testing for DirichletBC for mesh with 1-D line cells.
 */

#if !defined(pylith_bc_testdirichletbctet4_hh)
#define pylith_bc_testdirichletbcet4_hh

#include "TestDirichletBC.hh" // ISA TestDirichletBC

/// Namespace for pylith package
namespace pylith {
  namespace bc {
    class TestDirichletBCTet4;
  } // bc
} // pylith

/// C++ unit testing for DirichletBC for mesh with 3-D tet cells.
class pylith::bc::TestDirichletBCTet4 : public TestDirichletBC
{ // class TestDirichletBC

  // CPPUNIT TEST SUITE /////////////////////////////////////////////////
  CPPUNIT_TEST_SUB_SUITE( TestDirichletBCTet4, TestDirichletBC );

  CPPUNIT_TEST( testInitialize );
  CPPUNIT_TEST( testNumDimConstrained );
  CPPUNIT_TEST( testSetConstraintSizes );
  CPPUNIT_TEST( testSetConstraints );
  CPPUNIT_TEST( testSetField );
  CPPUNIT_TEST( testSetFieldIncr );

  CPPUNIT_TEST_SUITE_END();

  // PUBLIC METHODS /////////////////////////////////////////////////////
public :

  /// Setup testing data.
  void setUp(void);

}; // class TestDirichletBCTet4

#endif // pylith_bc_dirichletbctet4_hh


// End of file 
