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
 * @file unittests/libtests/bc/TestDirichletBCCases.hh
 *
 * @brief C++ TestDirichletBC object.
 *
 * Test cases C++ unit testing for DirichletBC for mesh.
 */

#if !defined(pylith_bc_testdirichletbccases_hh)
#define pylith_bc_testdirichletbccases_hh

#include "TestDirichletBC.hh" // ISA TestDirichletBC

/// Namespace for pylith package
namespace pylith {
  namespace bc {
    class TestDirichletBCTri3;
    class TestDirichletBCQuad4;
    class TestDirichletBCTet4;
    class TestDirichletBCHex8;
  } // bc
} // pylith

// ----------------------------------------------------------------------
/// C++ unit testing for DirichletBC for mesh with 2-D tri cells.
class pylith::bc::TestDirichletBCTri3 : public TestDirichletBC
{ // class TestDirichletBC

  // CPPUNIT TEST SUITE /////////////////////////////////////////////////
  CPPUNIT_TEST_SUB_SUITE( TestDirichletBCTri3, TestDirichletBC );

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

}; // class TestDirichletBCTri3


// ----------------------------------------------------------------------
/// C++ unit testing for DirichletBC for mesh with 2-D quad cells.
class pylith::bc::TestDirichletBCQuad4 : public TestDirichletBC
{ // class TestDirichletBC

  // CPPUNIT TEST SUITE /////////////////////////////////////////////////
  CPPUNIT_TEST_SUB_SUITE( TestDirichletBCQuad4, TestDirichletBC );

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

}; // class TestDirichletBCQuad4


// ----------------------------------------------------------------------
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


// ----------------------------------------------------------------------
/// C++ unit testing for DirichletBC for mesh with 3-D hex cells.
class pylith::bc::TestDirichletBCHex8 : public TestDirichletBC
{ // class TestDirichletBC

  // CPPUNIT TEST SUITE /////////////////////////////////////////////////
  CPPUNIT_TEST_SUB_SUITE( TestDirichletBCHex8, TestDirichletBC );

  CPPUNIT_TEST( testInitialize );
  CPPUNIT_TEST( testSetConstraintSizes );
  CPPUNIT_TEST( testSetConstraints );
  CPPUNIT_TEST( testSetField );
  CPPUNIT_TEST( testSetFieldIncr );

  CPPUNIT_TEST_SUITE_END();

  // PUBLIC METHODS /////////////////////////////////////////////////////
public :

  /// Setup testing data.
  void setUp(void);

}; // class TestDirichletBCHex8


#endif // pylith_bc_dirichletbccases_hh


// End of file 
