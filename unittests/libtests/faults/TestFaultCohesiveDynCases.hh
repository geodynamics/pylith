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
 * @file unittests/libtests/faults/TestFaultCohesiveDynCases.hh
 *
 * @brief C++ TestFaultCohesiveDynCaes object.
 *
 * C++ unit testing for FaultCohesiveDyn.
 */

#if !defined(pylith_faults_testfaultcohesivedyncases_hh)
#define pylith_faults_testfaultcohesivedyncases_hh

#include "TestFaultCohesiveDyn.hh" // ISA TestFaultCohesiveDyn

/// Namespace for pylith package
namespace pylith {
  namespace faults {
    class TestFaultCohesiveDynTri3;
    class TestFaultCohesiveDynTri3d;
    class TestFaultCohesiveDynQuad4;
    class TestFaultCohesiveDynTet4;
    class TestFaultCohesiveDynHex8;
  } // bc
} // pylith

// ----------------------------------------------------------------------
/// C++ unit testing for FaultCohesiveDyn for mesh with 2-D triangular cells.
class pylith::faults::TestFaultCohesiveDynTri3 : public TestFaultCohesiveDyn
{ // class TestFaultCohesiveDynTri3

  // CPPUNIT TEST SUITE /////////////////////////////////////////////////
  CPPUNIT_TEST_SUITE( TestFaultCohesiveDynTri3 );

  CPPUNIT_TEST( testInitialize );
  CPPUNIT_TEST( testConstrainSolnSpaceStick );
  CPPUNIT_TEST( testConstrainSolnSpaceSlip );
  CPPUNIT_TEST( testConstrainSolnSpaceOpen );
  CPPUNIT_TEST( testUpdateStateVars );
  CPPUNIT_TEST( testCalcTractions );

  CPPUNIT_TEST_SUITE_END();

  // PUBLIC METHODS /////////////////////////////////////////////////////
public :

  /// Setup testing data.
  void setUp(void);

}; // class TestFaultCohesiveDynTri3


// ----------------------------------------------------------------------
/// C++ unit testing for FaultCohesiveDyn for mesh with 2-D triangular cells.
class pylith::faults::TestFaultCohesiveDynTri3d : public TestFaultCohesiveDyn
{ // class TestFaultCohesiveDynTri3d

  // CPPUNIT TEST SUITE /////////////////////////////////////////////////
  CPPUNIT_TEST_SUITE( TestFaultCohesiveDynTri3d );

  CPPUNIT_TEST( testInitialize );
  CPPUNIT_TEST( testConstrainSolnSpaceStick );
  CPPUNIT_TEST( testConstrainSolnSpaceSlip );
  CPPUNIT_TEST( testConstrainSolnSpaceOpen );
  CPPUNIT_TEST( testUpdateStateVars );
  CPPUNIT_TEST( testCalcTractions );

  CPPUNIT_TEST_SUITE_END();

  // PUBLIC METHODS /////////////////////////////////////////////////////
public :

  /// Setup testing data.
  void setUp(void);

}; // class TestFaultCohesiveDynTri3d


// ----------------------------------------------------------------------
/// C++ unit testing for FaultCohesiveDyn for mesh with 2-D quadrilateral cells.
class pylith::faults::TestFaultCohesiveDynQuad4 : public TestFaultCohesiveDyn
{ // class TestFaultCohesiveDynQuad4

  // CPPUNIT TEST SUITE /////////////////////////////////////////////////
  CPPUNIT_TEST_SUITE( TestFaultCohesiveDynQuad4 );

  CPPUNIT_TEST( testInitialize );
  CPPUNIT_TEST( testConstrainSolnSpaceStick );
  CPPUNIT_TEST( testConstrainSolnSpaceSlip );
  CPPUNIT_TEST( testConstrainSolnSpaceOpen );
  CPPUNIT_TEST( testUpdateStateVars );
  CPPUNIT_TEST( testCalcTractions );

  CPPUNIT_TEST_SUITE_END();

  // PUBLIC METHODS /////////////////////////////////////////////////////
public :

  /// Setup testing data.
  void setUp(void);

}; // class TestFaultCohesiveDynQuad4


// ----------------------------------------------------------------------
/// C++ unit testing for FaultCohesiveDyn for mesh with 2-D quadrilateral cells.
class pylith::faults::TestFaultCohesiveDynTet4 : public TestFaultCohesiveDyn
{ // class TestFaultCohesiveDynTet4

  // CPPUNIT TEST SUITE /////////////////////////////////////////////////
  CPPUNIT_TEST_SUITE( TestFaultCohesiveDynTet4 );

  CPPUNIT_TEST( testInitialize );
  CPPUNIT_TEST( testConstrainSolnSpaceStick );
  CPPUNIT_TEST( testConstrainSolnSpaceSlip );
  CPPUNIT_TEST( testConstrainSolnSpaceOpen );
  CPPUNIT_TEST( testUpdateStateVars );
  CPPUNIT_TEST( testCalcTractions );

  CPPUNIT_TEST_SUITE_END();

  // PUBLIC METHODS /////////////////////////////////////////////////////
public :

  /// Setup testing data.
  void setUp(void);

}; // class TestFaultCohesiveDynTet4


// ----------------------------------------------------------------------
/// C++ unit testing for FaultCohesiveDyn for mesh with 3-D quadrilateral cells.
class pylith::faults::TestFaultCohesiveDynHex8 : public TestFaultCohesiveDyn
{ // class TestFaultCohesiveDynHex8

  // CPPUNIT TEST SUITE /////////////////////////////////////////////////
  CPPUNIT_TEST_SUITE( TestFaultCohesiveDynHex8 );

  CPPUNIT_TEST( testInitialize );
  CPPUNIT_TEST( testConstrainSolnSpaceStick );
  CPPUNIT_TEST( testConstrainSolnSpaceSlip );
  CPPUNIT_TEST( testConstrainSolnSpaceOpen );
  CPPUNIT_TEST( testUpdateStateVars );
  CPPUNIT_TEST( testCalcTractions );

  CPPUNIT_TEST_SUITE_END();

  // PUBLIC METHODS /////////////////////////////////////////////////////
public :

  /// Setup testing data.
  void setUp(void);

}; // class TestFaultCohesiveDynHex8


#endif // pylith_faults_testfaultcohesivedyncases_hh


// End of file 
