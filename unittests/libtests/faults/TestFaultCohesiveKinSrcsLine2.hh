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
 * @file unittests/libtests/faults/TestFaultCohesiveKinSrcsLine2.hh
 *
 * @brief C++ TestFaultCohesiveKinSrcsLine2 object.
 *
 * C++ unit testing for FaultCohesiveKin for mesh with 1-D line cells
 * and multiple earthquake ruptures.
 */

#if !defined(pylith_faults_testfaultcohesivekinsrcsline2_hh)
#define pylith_faults_testfaultcohesivekinsrcsline2_hh

#include "TestFaultCohesiveKinSrcs.hh" // ISA TestFaultCohesiveKinSrcs

/// Namespace for pylith package
namespace pylith {
  namespace faults {
    class TestFaultCohesiveKinSrcsLine2;
  } // bc
} // pylith

/// C++ unit testing for FaultCohesiveKin for mesh with 1-D line cells.
class pylith::faults::TestFaultCohesiveKinSrcsLine2 : public TestFaultCohesiveKinSrcs
{ // class TestFaultCohesiveKinSrcsLine2

  // CPPUNIT TEST SUITE /////////////////////////////////////////////////
  CPPUNIT_TEST_SUITE( TestFaultCohesiveKinSrcsLine2 );

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

}; // class TestFaultCohesiveKinSrcsLine2

#endif // pylith_faults_testfaultcohesivesrcsline2_hh


// End of file 
