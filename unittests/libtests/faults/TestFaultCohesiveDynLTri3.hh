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
 * @file unittests/libtests/faults/TestFaultCohesiveDynLTri3.hh
 *
 * @brief C++ TestFaultCohesiveDynLTri3 object.
 *
 * C++ unit testing for FaultCohesiveDynL for mesh with 2-D triangular cells.
 */

#if !defined(pylith_faults_testfaultcohesivedynltri3_hh)
#define pylith_faults_testfaultcohesivedynltri3_hh

#include "TestFaultCohesiveDynL.hh" // ISA TestFaultCohesiveDynL

/// Namespace for pylith package
namespace pylith {
  namespace faults {
    class TestFaultCohesiveDynLTri3;
  } // bc
} // pylith

/// C++ unit testing for FaultCohesiveDynL for mesh with 2-D triangular cells.
class pylith::faults::TestFaultCohesiveDynLTri3 : public TestFaultCohesiveDynL
{ // class TestFaultCohesiveDynLTri3

  // CPPUNIT TEST SUITE /////////////////////////////////////////////////
  CPPUNIT_TEST_SUITE( TestFaultCohesiveDynLTri3 );

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

}; // class TestFaultCohesiveDynLTri3

#endif // pylith_faults_testfaultcohesivedynltri3_hh


// End of file 
