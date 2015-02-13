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
 * @file unittests/libtests/bc/TestDirichletBCMultiCases.hh
 *
 * @brief C++ TestDirichletBC object.
 *
 * Test cases for C++ unit testing for DirichletBC for mesh.
 */

#if !defined(pylith_bc_testdirichletbcmulticases_hh)
#define pylith_bc_testdirichletbcmulticases_hh

#include "TestDirichletBCMulti.hh" // ISA TestDirichletBC

/// Namespace for pylith package
namespace pylith {
  namespace bc {
    class TestDirichletBCMultiTri3;
    class TestDirichletBCMultiTet4;
  } // bc
} // pylith

// ----------------------------------------------------------------------
/// C++ unit testing for DirichletBC for mesh with 2-D tri cells.
class pylith::bc::TestDirichletBCMultiTri3 : public TestDirichletBCMulti
{ // class TestDirichletBC

  // CPPUNIT TEST SUITE /////////////////////////////////////////////////
  CPPUNIT_TEST_SUITE( TestDirichletBCMultiTri3 );
  CPPUNIT_TEST( testSetConstraintSizes );
  CPPUNIT_TEST( testSetConstraints );
  CPPUNIT_TEST( testSetField );
  CPPUNIT_TEST_SUITE_END();

  // PUBLIC METHODS /////////////////////////////////////////////////////
public :

  /// Setup testing data.
  void setUp(void);

}; // class TestDirichletBCMultiTri3


// ----------------------------------------------------------------------
/// C++ unit testing for DirichletBC for mesh with 3-D tet cells.
class pylith::bc::TestDirichletBCMultiTet4 : public TestDirichletBCMulti
{ // class TestDirichletBC

  // CPPUNIT TEST SUITE /////////////////////////////////////////////////
  CPPUNIT_TEST_SUITE( TestDirichletBCMultiTet4 );

  CPPUNIT_TEST( testSetConstraintSizes );
  CPPUNIT_TEST( testSetConstraints );
  CPPUNIT_TEST( testSetField );
  CPPUNIT_TEST( testSetFieldIncr );

  CPPUNIT_TEST_SUITE_END();

  // PUBLIC METHODS /////////////////////////////////////////////////////
public :

  /// Setup testing data.
  void setUp(void);

}; // class TestDirichletBCMultiTet4


#endif // pylith_bc_dirichletbcmulticases_hh


// End of file 
