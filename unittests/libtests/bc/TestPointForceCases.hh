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
 * @file unittests/libtests/bc/TestPointForceCases.hh
 *
 * @brief C++ TestPointForce object.
 *
 * Test cases for C++ unit testing for PointForce.
 */

#if !defined(pylith_bc_testpointforcebccases_hh)
#define pylith_bc_testpointforcebccases_hh

#include "TestPointForce.hh" // ISA TestPointForce

/// Namespace for pylith package
namespace pylith {
  namespace bc {
    class TestPointForceTri3;
    class TestPointForceQuad4;
    class TestPointForceTet4;
    class TestPointForceHex8;
  } // bc
} // pylith

// ----------------------------------------------------------------------
/// C++ unit testing for PointForce for mesh with 2-D tri cells.
class pylith::bc::TestPointForceTri3 : public TestPointForce
{ // class TestPointForce

  // CPPUNIT TEST SUITE /////////////////////////////////////////////////
  CPPUNIT_TEST_SUB_SUITE( TestPointForceTri3, TestPointForce );

  CPPUNIT_TEST( testInitialize );
  CPPUNIT_TEST( testIntegrateResidual );
  CPPUNIT_TEST( testVerifyConfiguration );

  CPPUNIT_TEST_SUITE_END();

  // PUBLIC METHODS /////////////////////////////////////////////////////
public :

  /// Setup testing data.
  void setUp(void);

}; // class TestPointForceTri3


// ----------------------------------------------------------------------
/// C++ unit testing for PointForce for mesh with 2-D quad4 cells.
class pylith::bc::TestPointForceQuad4 : public TestPointForce
{ // class TestPointForce

  // CPPUNIT TEST SUITE /////////////////////////////////////////////////
  CPPUNIT_TEST_SUB_SUITE( TestPointForceQuad4, TestPointForce );

  CPPUNIT_TEST( testInitialize );
  CPPUNIT_TEST( testIntegrateResidual );
  CPPUNIT_TEST( testVerifyConfiguration );

  CPPUNIT_TEST_SUITE_END();

  // PUBLIC METHODS /////////////////////////////////////////////////////
public :

  /// Setup testing data.
  void setUp(void);

}; // class TestPointForceQuad4

// ----------------------------------------------------------------------
/// C++ unit testing for PointForce for mesh with 3-D tet cells.
class pylith::bc::TestPointForceTet4 : public TestPointForce
{ // class TestPointForce

  // CPPUNIT TEST SUITE /////////////////////////////////////////////////
  CPPUNIT_TEST_SUB_SUITE( TestPointForceTet4, TestPointForce );

  CPPUNIT_TEST( testInitialize );
  CPPUNIT_TEST( testIntegrateResidual );
  CPPUNIT_TEST( testVerifyConfiguration );

  CPPUNIT_TEST_SUITE_END();

  // PUBLIC METHODS /////////////////////////////////////////////////////
public :

  /// Setup testing data.
  void setUp(void);

}; // class TestPointForceTet4


// ----------------------------------------------------------------------
/// C++ unit testing for PointForce for mesh with 3-D hex cells.
class pylith::bc::TestPointForceHex8 : public TestPointForce
{ // class TestPointForce

  // CPPUNIT TEST SUITE /////////////////////////////////////////////////
  CPPUNIT_TEST_SUB_SUITE( TestPointForceHex8, TestPointForce );

  CPPUNIT_TEST( testInitialize );
  CPPUNIT_TEST( testIntegrateResidual );
  CPPUNIT_TEST( testVerifyConfiguration );

  CPPUNIT_TEST_SUITE_END();

  // PUBLIC METHODS /////////////////////////////////////////////////////
public :

  /// Setup testing data.
  void setUp(void);

}; // class TestPointForceHex8


#endif // pylith_bc_pointforcebcline2_hh


// End of file 
