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
// Copyright (c) 2010-2011 University of California, Davis
//
// See COPYING for license information.
//
// ----------------------------------------------------------------------
//

/**
 * @file unittests/libtests/faults/TestFaultCohesiveKinTri3.hh
 *
 * @brief C++ TestFaultCohesiveKinTri3 object.
 *
 * C++ unit testing for FaultCohesiveKin for mesh with 2-D triangular cells.
 */

#if !defined(pylith_faults_testfaultcohesivekintri3_hh)
#define pylith_faults_testfaultcohesivekintri3_hh

#include "TestFaultCohesiveKin.hh" // ISA TestFaultCohesiveKin

/// Namespace for pylith package
namespace pylith {
  namespace faults {
    class TestFaultCohesiveKinTri3;
  } // bc
} // pylith

/// C++ unit testing for FaultCohesiveKin for mesh with 2-D triangular cells.
class pylith::faults::TestFaultCohesiveKinTri3 : public TestFaultCohesiveKin
{ // class TestFaultCohesiveKinTri3

  // CPPUNIT TEST SUITE /////////////////////////////////////////////////
  CPPUNIT_TEST_SUITE( TestFaultCohesiveKinTri3 );

  CPPUNIT_TEST( testInitialize );
  CPPUNIT_TEST( testIntegrateResidual );
  CPPUNIT_TEST( testIntegrateJacobian );
  CPPUNIT_TEST( testIntegrateJacobianLumped );
  CPPUNIT_TEST( testAdjustSolnLumped );
  CPPUNIT_TEST( testUpdateStateVars );
  CPPUNIT_TEST( testCalcTractionsChange );
  CPPUNIT_TEST( testSplitField );

  CPPUNIT_TEST_SUITE_END();

  // PUBLIC METHODS /////////////////////////////////////////////////////
public :

  /// Setup testing data.
  void setUp(void);

}; // class TestFaultCohesiveKinTri3

#endif // pylith_faults_testfaultcohesivetri3_hh


// End of file 
