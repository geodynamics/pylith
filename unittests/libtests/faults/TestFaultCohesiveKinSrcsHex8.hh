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
// Copyright (c) 2010-2013 University of California, Davis
//
// See COPYING for license information.
//
// ----------------------------------------------------------------------
//

/**
 * @file unittests/libtests/faults/TestFaultCohesiveKinSrcsHex8.hh
 *
 * @brief C++ TestFaultCohesiveKinSrcsHex8 object.
 *
 * C++ unit testing for FaultCohesiveKin for mesh with 3-D hex cells.
 */

#if !defined(pylith_faults_testfaultcohesivekinsrcshex8_hh)
#define pylith_faults_testfaultcohesivekinsrcshex8_hh

#include "TestFaultCohesiveKinSrcs.hh" // ISA TestFaultCohesiveKinSrcs

/// Namespace for pylith package
namespace pylith {
  namespace faults {
    class TestFaultCohesiveKinSrcsHex8;
  } // bc
} // pylith

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

#endif // pylith_faults_testfaultcohesivekinsrcshex8_hh


// End of file 
