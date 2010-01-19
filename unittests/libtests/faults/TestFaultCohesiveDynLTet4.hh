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
 * @file unittests/libtests/faults/TestFaultCohesiveDynLTet4.hh
 *
 * @brief C++ TestFaultCohesiveDynLTet4 object.
 *
 * C++ unit testing for FaultCohesiveDynL for mesh with 2-D quadrilateral cells.
 */

#if !defined(pylith_faults_testfaultcohesivedynltet4_hh)
#define pylith_faults_testfaultcohesivedynltet4_hh

#include "TestFaultCohesiveDynL.hh" // ISA TestFaultCohesiveDynL

/// Namespace for pylith package
namespace pylith {
  namespace faults {
    class TestFaultCohesiveDynLTet4;
  } // bc
} // pylith

/// C++ unit testing for FaultCohesiveDynL for mesh with 2-D quadrilateral cells.
class pylith::faults::TestFaultCohesiveDynLTet4 : public TestFaultCohesiveDynL
{ // class TestFaultCohesiveDynLTet4

  // CPPUNIT TEST SUITE /////////////////////////////////////////////////
  CPPUNIT_TEST_SUITE( TestFaultCohesiveDynLTet4 );

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

}; // class TestFaultCohesiveDynLTet4

#endif // pylith_faults_testfaultcohesivedynltet4_hh


// End of file 
