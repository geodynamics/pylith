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
 * @file unittests/libtests/faults/TestFaultCohesiveDynTet4.hh
 *
 * @brief C++ TestFaultCohesiveDynTet4 object.
 *
 * C++ unit testing for FaultCohesiveDyn for mesh with 2-D quadrilateral cells.
 */

#if !defined(pylith_faults_testfaultcohesivedyntet4_hh)
#define pylith_faults_testfaultcohesivedyntet4_hh

#include "TestFaultCohesiveDyn.hh" // ISA TestFaultCohesiveDyn

/// Namespace for pylith package
namespace pylith {
  namespace faults {
    class TestFaultCohesiveDynTet4;
  } // bc
} // pylith

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

#endif // pylith_faults_testfaultcohesivedyntet4_hh


// End of file 
