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
 * @file unittests/libtests/faults/TestFaultCohesiveDynTri3.hh
 *
 * @brief C++ TestFaultCohesiveDynTri3 object.
 *
 * C++ unit testing for FaultCohesiveDyn for mesh with 2-D triangular cells.
 */

#if !defined(pylith_faults_testfaultcohesivedyntri3_hh)
#define pylith_faults_testfaultcohesivedyntri3_hh

#include "TestFaultCohesiveDyn.hh" // ISA TestFaultCohesiveDyn

/// Namespace for pylith package
namespace pylith {
  namespace faults {
    class TestFaultCohesiveDynTri3;
  } // bc
} // pylith

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

#endif // pylith_faults_testfaultcohesivedyntri3_hh


// End of file 
