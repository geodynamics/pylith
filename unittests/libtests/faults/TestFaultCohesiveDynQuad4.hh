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
 * @file unittests/libtests/faults/TestFaultCohesiveDynQuad4.hh
 *
 * @brief C++ TestFaultCohesiveDynQuad4 object.
 *
 * C++ unit testing for FaultCohesiveDyn for mesh with 2-D quadrilateral cells.
 */

#if !defined(pylith_faults_testfaultcohesivedynquad4_hh)
#define pylith_faults_testfaultcohesivedynquad4_hh

#include "TestFaultCohesiveDyn.hh" // ISA TestFaultCohesiveDyn

/// Namespace for pylith package
namespace pylith {
  namespace faults {
    class TestFaultCohesiveDynQuad4;
  } // bc
} // pylith

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

#endif // pylith_faults_testfaultcohesivedynquad4_hh


// End of file 
