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
 * @file unittests/libtests/bc/TestAbsorbingDampersCases.hh
 *
 * @brief C++ TestAbsorbingDampers object.
 *
 * Test cases for C++ unit testing for AbsorbingDampers.
 */

#if !defined(pylith_bc_testabsorbingdamperscases_hh)
#define pylith_bc_testabsorbingdamperscases_hh

#include "TestAbsorbingDampers.hh" // ISA TestAbsorbingDampers

/// Namespace for pylith package
namespace pylith {
  namespace bc {
    class TestAbsorbingDampersTri3;
    class TestAbsorbingDampersQuad4;
    class TestAbsorbingDampersTet4;
    class TestAbsorbingDampersHex8;
  } // bc
} // pylith

// ----------------------------------------------------------------------
/// C++ unit testing for AbsorbingDampers for mesh with 2-D tri cells.
class pylith::bc::TestAbsorbingDampersTri3 : public TestAbsorbingDampers
{ // class TestAbsorbingDampers

  // CPPUNIT TEST SUITE /////////////////////////////////////////////////
  CPPUNIT_TEST_SUB_SUITE( TestAbsorbingDampersTri3, TestAbsorbingDampers );
  CPPUNIT_TEST( testInitialize );
  CPPUNIT_TEST( testIntegrateResidual );
  CPPUNIT_TEST( testIntegrateJacobian );
  CPPUNIT_TEST( testIntegrateJacobianLumped );
  CPPUNIT_TEST_SUITE_END();

  // PUBLIC METHODS /////////////////////////////////////////////////////
public :

  /// Setup testing data.
  void setUp(void);

}; // class TestAbsorbingDampersTri3

// ----------------------------------------------------------------------
/// C++ unit testing for AbsorbingDampers for mesh with 2-D quad cells.
class pylith::bc::TestAbsorbingDampersQuad4 : public TestAbsorbingDampers
{ // class TestAbsorbingDampers

  // CPPUNIT TEST SUITE /////////////////////////////////////////////////
  CPPUNIT_TEST_SUB_SUITE( TestAbsorbingDampersQuad4, TestAbsorbingDampers );
  CPPUNIT_TEST( testInitialize );
  CPPUNIT_TEST( testIntegrateResidual );
  CPPUNIT_TEST( testIntegrateJacobian );
  CPPUNIT_TEST( testIntegrateJacobianLumped );
  CPPUNIT_TEST_SUITE_END();

  // PUBLIC METHODS /////////////////////////////////////////////////////
public :

  /// Setup testing data.
  void setUp(void);

}; // class TestAbsorbingDampersQuad4


// ----------------------------------------------------------------------
/// C++ unit testing for AbsorbingDampers for mesh with 3-D tet cells.
class pylith::bc::TestAbsorbingDampersTet4 : public TestAbsorbingDampers
{ // class TestAbsorbingDampers

  // CPPUNIT TEST SUITE /////////////////////////////////////////////////
  CPPUNIT_TEST_SUB_SUITE( TestAbsorbingDampersTet4, TestAbsorbingDampers );
  CPPUNIT_TEST( testInitialize );
  CPPUNIT_TEST( testIntegrateResidual );
  CPPUNIT_TEST( testIntegrateJacobian );
  CPPUNIT_TEST( testIntegrateJacobianLumped );
  CPPUNIT_TEST_SUITE_END();

  // PUBLIC METHODS /////////////////////////////////////////////////////
public :

  /// Setup testing data.
  void setUp(void);

}; // class TestAbsorbingDampersTet4

// ----------------------------------------------------------------------
/// C++ unit testing for AbsorbingDampers for mesh with 3-D hex cells.
class pylith::bc::TestAbsorbingDampersHex8 : public TestAbsorbingDampers
{ // class TestAbsorbingDampers

  // CPPUNIT TEST SUITE /////////////////////////////////////////////////
  CPPUNIT_TEST_SUB_SUITE( TestAbsorbingDampersHex8, TestAbsorbingDampers );

  CPPUNIT_TEST( testInitialize );
  CPPUNIT_TEST( testIntegrateResidual );
  CPPUNIT_TEST( testIntegrateJacobian );
  CPPUNIT_TEST( testIntegrateJacobianLumped );

  CPPUNIT_TEST_SUITE_END();

  // PUBLIC METHODS /////////////////////////////////////////////////////
public :

  /// Setup testing data.
  void setUp(void);

}; // class TestAbsorbingDampersHex8

#endif // pylith_bc_absorbingdamperscases_hh


// End of file 
