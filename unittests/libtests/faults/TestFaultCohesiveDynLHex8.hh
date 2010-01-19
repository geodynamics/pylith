// -*- C++ -*-
//
// ----------------------------------------------------------------------
//
//                           Brad T. Aagaard
//                        U.S. Geological Survey
//
// {LicenseText}
//
// ----------------------------------------------------------------------
//

/**
 * @file unittests/libtests/faults/TestFaultCohesiveDynLHex8.hh
 *
 * @brief C++ TestFaultCohesiveDynLHex8 object.
 *
 * C++ unit testing for FaultCohesiveDynL for mesh with 3-D quadrilateral cells.
 */

#if !defined(pylith_faults_testfaultcohesivedynlhex8_hh)
#define pylith_faults_testfaultcohesivedynlhex8_hh

#include "TestFaultCohesiveDynL.hh" // ISA TestFaultCohesiveDynL

/// Namespace for pylith package
namespace pylith {
  namespace faults {
    class TestFaultCohesiveDynLHex8;
  } // bc
} // pylith

/// C++ unit testing for FaultCohesiveDynL for mesh with 3-D quadrilateral cells.
class pylith::faults::TestFaultCohesiveDynLHex8 : public TestFaultCohesiveDynL
{ // class TestFaultCohesiveDynLHex8

  // CPPUNIT TEST SUITE /////////////////////////////////////////////////
  CPPUNIT_TEST_SUITE( TestFaultCohesiveDynLHex8 );

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

}; // class TestFaultCohesiveDynLHex8

#endif // pylith_faults_testfaultcohesivedynlhex8_hh


// End of file 
