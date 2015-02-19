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
 * @file unittests/libtests/faults/TestFaultCohesiveKinSrcsCases.hh
 *
 * @brief C++ unit testing for FaultCohesiveKin with multiple earthquake ruptures.
 */

#if !defined(pylith_faults_testfaultcohesivekinsrcscases_hh)
#define pylith_faults_testfaultcohesivekinsrcscases_hh

#include "TestFaultCohesiveKinSrcs.hh" // ISA TestFaultCohesiveKinSrcs

/// Namespace for pylith package
namespace pylith {
  namespace faults {
    class TestFaultCohesiveKinSrcsTri3;
    class TestFaultCohesiveKinSrcsQuad4;
    class TestFaultCohesiveKinSrcsTet4;
    class TestFaultCohesiveKinSrcsHex8;
  } // bc
} // pylith

// ----------------------------------------------------------------------
/// C++ unit testing for FaultCohesiveKinSrcs for mesh with 2-D triangular cells.
class pylith::faults::TestFaultCohesiveKinSrcsTri3 : public TestFaultCohesiveKinSrcs
{ // class TestFaultCohesiveKinSrcsTri3

  // CPPUNIT TEST SUITE /////////////////////////////////////////////////
  CPPUNIT_TEST_SUITE( TestFaultCohesiveKinSrcsTri3 );

  CPPUNIT_TEST( testInitialize );
  CPPUNIT_TEST( testIntegrateResidual );
  CPPUNIT_TEST( testIntegrateJacobian );
  CPPUNIT_TEST( testIntegrateJacobianLumped );
  CPPUNIT_TEST( testCalcTractionsChange );

  CPPUNIT_TEST_SUITE_END();

  // PUBLIC METHODS /////////////////////////////////////////////////////
public :

  /// Setup testing data.
  void setUp(void);

}; // class TestFaultCohesiveKinSrcsTri3


// ----------------------------------------------------------------------
/// C++ unit testing for FaultCohesiveKin for mesh with 2-D quadrilateral cells.
class pylith::faults::TestFaultCohesiveKinSrcsQuad4 : public TestFaultCohesiveKinSrcs
{ // class TestFaultCohesiveKinSrcsQuad4

  // CPPUNIT TEST SUITE /////////////////////////////////////////////////
  CPPUNIT_TEST_SUITE( TestFaultCohesiveKinSrcsQuad4 );

  CPPUNIT_TEST( testInitialize );
  CPPUNIT_TEST( testIntegrateResidual );
  CPPUNIT_TEST( testIntegrateJacobian );
  CPPUNIT_TEST( testIntegrateJacobianLumped );
  CPPUNIT_TEST( testCalcTractionsChange );

  CPPUNIT_TEST_SUITE_END();

  // PUBLIC METHODS /////////////////////////////////////////////////////
public :

  /// Setup testing data.
  void setUp(void);

}; // class TestFaultCohesiveKinSrcsQuad4


// ----------------------------------------------------------------------
/// C++ unit testing for FaultCohesiveKin for mesh with 3-D tetrahedral cells.
class pylith::faults::TestFaultCohesiveKinSrcsTet4 : public TestFaultCohesiveKinSrcs
{ // class TestFaultCohesiveKinSrcsTet4

  // CPPUNIT TEST SUITE /////////////////////////////////////////////////
  CPPUNIT_TEST_SUITE( TestFaultCohesiveKinSrcsTet4 );

  CPPUNIT_TEST( testInitialize );
  CPPUNIT_TEST( testIntegrateResidual );
  CPPUNIT_TEST( testIntegrateJacobian );
  CPPUNIT_TEST( testIntegrateJacobianLumped );
  CPPUNIT_TEST( testCalcTractionsChange );

  CPPUNIT_TEST_SUITE_END();

  // PUBLIC METHODS /////////////////////////////////////////////////////
public :

  /// Setup testing data.
  void setUp(void);

}; // class TestFaultCohesiveKinSrcsTet4


// ----------------------------------------------------------------------
/// C++ unit testing for FaultCohesiveKin for mesh with 3-D hex cells.
class pylith::faults::TestFaultCohesiveKinSrcsHex8 : public TestFaultCohesiveKinSrcs
{ // class TestFaultCohesiveKinSrcsHex8

  // CPPUNIT TEST SUITE /////////////////////////////////////////////////
  CPPUNIT_TEST_SUITE( TestFaultCohesiveKinSrcsHex8 );

  CPPUNIT_TEST( testInitialize );
  CPPUNIT_TEST( testIntegrateResidual );
  CPPUNIT_TEST( testIntegrateJacobian );
  CPPUNIT_TEST( testIntegrateJacobianLumped );
  CPPUNIT_TEST( testCalcTractionsChange );

  CPPUNIT_TEST_SUITE_END();

  // PUBLIC METHODS /////////////////////////////////////////////////////
public :

  /// Setup testing data.
  void setUp(void);

}; // class TestFaultCohesiveKinSrcsHex8


#endif // pylith_faults_testfaultcohesivesrcscases_hh


// End of file 
