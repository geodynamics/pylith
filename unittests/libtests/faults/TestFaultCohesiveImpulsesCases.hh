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
 * @file unittests/libtests/faults/TestFaultCohesiveImpulsesCases.hh
 *
 * @brief C++ unit testing for FaultCohesiveImpulses with various cell
 * types.
 */

#if !defined(pylith_faults_testfaultcohesiveimpulsescases_hh)
#define pylith_faults_testfaultcohesiveimpulsescases_hh

#include "TestFaultCohesiveImpulses.hh"

/// Namespace for pylith package
namespace pylith {
  namespace faults {
    class TestFaultCohesiveImpulsesTri3;
    class TestFaultCohesiveImpulsesQuad4;
    class TestFaultCohesiveImpulsesTet4;
    class TestFaultCohesiveImpulsesHex8;
  } // faults
} // pylith


// ----------------------------------------------------------------------
/// C++ unit testing for FaultCohesiveImpulses with tri3 cells
class pylith::faults::TestFaultCohesiveImpulsesTri3 : public TestFaultCohesiveImpulses
{ // class TestFaultCohesiveImpulsesTri3

  // CPPUNIT TEST SUITE /////////////////////////////////////////////////
  CPPUNIT_TEST_SUITE( TestFaultCohesiveImpulsesTri3 );

  CPPUNIT_TEST( testNumImpulses );
  CPPUNIT_TEST( testInitialize );
  CPPUNIT_TEST( testIntegrateResidual );

  CPPUNIT_TEST_SUITE_END();

  // PUBLIC METHODS /////////////////////////////////////////////////////
public :

  /// Setup testing data.
  void setUp(void);

}; // class TestFaultCohesiveImpulsesTri3


// ----------------------------------------------------------------------
/// C++ unit testing for FaultCohesiveImpulses with Quad4 cells
class pylith::faults::TestFaultCohesiveImpulsesQuad4 : public TestFaultCohesiveImpulses
{ // class TestFaultCohesiveImpulsesQuad4

  // CPPUNIT TEST SUITE /////////////////////////////////////////////////
  CPPUNIT_TEST_SUITE( TestFaultCohesiveImpulsesQuad4 );

  CPPUNIT_TEST( testNumImpulses );
  CPPUNIT_TEST( testInitialize );
  CPPUNIT_TEST( testIntegrateResidual );

  CPPUNIT_TEST_SUITE_END();

  // PUBLIC METHODS /////////////////////////////////////////////////////
public :

  /// Setup testing data.
  void setUp(void);

}; // class TestFaultCohesiveImpulsesQuad4


// ----------------------------------------------------------------------
/// C++ unit testing for FaultCohesiveImpulses with Tet4 cells
class pylith::faults::TestFaultCohesiveImpulsesTet4 : public TestFaultCohesiveImpulses
{ // class TestFaultCohesiveImpulsesTet4

  // CPPUNIT TEST SUITE /////////////////////////////////////////////////
  CPPUNIT_TEST_SUITE( TestFaultCohesiveImpulsesTet4 );

  CPPUNIT_TEST( testNumImpulses );
  CPPUNIT_TEST( testInitialize );
  CPPUNIT_TEST( testIntegrateResidual );

  CPPUNIT_TEST_SUITE_END();

  // PUBLIC METHODS /////////////////////////////////////////////////////
public :

  /// Setup testing data.
  void setUp(void);

}; // class TestFaultCohesiveImpulsesTet4


// ----------------------------------------------------------------------
/// C++ unit testing for FaultCohesiveImpulses with Hex8 cells
class pylith::faults::TestFaultCohesiveImpulsesHex8 : public TestFaultCohesiveImpulses
{ // class TestFaultCohesiveImpulsesHex8

  // CPPUNIT TEST SUITE /////////////////////////////////////////////////
  CPPUNIT_TEST_SUITE( TestFaultCohesiveImpulsesHex8 );

  CPPUNIT_TEST( testNumImpulses );
  CPPUNIT_TEST( testInitialize );
  CPPUNIT_TEST( testIntegrateResidual );

  CPPUNIT_TEST_SUITE_END();

  // PUBLIC METHODS /////////////////////////////////////////////////////
public :

  /// Setup testing data.
  void setUp(void);

}; // class TestFaultCohesiveImpulsesHex8



#endif // pylith_faults_testfaultcohesiveimpulsescases_hh


// End of file 
