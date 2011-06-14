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
 * @file unittests/libtests/faults/TestFaultCohesiveDynTri3d.hh
 *
 * @brief C++ TestFaultCohesiveDynTri3d object.
 *
 * C++ unit testing for FaultCohesiveDyn for mesh with 2-D triangular cells.
 */

#if !defined(pylith_faults_testfaultcohesivedyntri3d_hh)
#define pylith_faults_testfaultcohesivedyntri3d_hh

#include "TestFaultCohesiveDyn.hh" // ISA TestFaultCohesiveDyn

/// Namespace for pylith package
namespace pylith {
  namespace faults {
    class TestFaultCohesiveDynTri3d;
  } // bc
} // pylith

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

#endif // pylith_faults_testfaultcohesivedyntri3d_hh


// End of file 
