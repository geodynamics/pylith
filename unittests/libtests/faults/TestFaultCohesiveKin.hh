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
 * @file unittests/libtests/faults/TestFaultCohesiveKin.hh
 *
 * @brief C++ TestFaultCohesiveKin object
 *
 * C++ unit testing for FaultCohesiveKin.
 */

#if !defined(pylith_faults_testfaultcohesivekin_hh)
#define pylith_faults_testfaultcohesivekin_hh

#include "TestFaultCohesive.hh"

/// Namespace for pylith package
namespace pylith {
  namespace faults {
    class TestFaultCohesiveKin;
  } // faults
} // pylith

/// C++ unit testing for FaultCohesiveKin
class pylith::faults::TestFaultCohesiveKin : public TestFaultCohesive
{ // class TestFaultCohesiveKin

  // CPPUNIT TEST SUITE /////////////////////////////////////////////////
  CPPUNIT_TEST_SUITE( TestFaultCohesiveKin );
  CPPUNIT_TEST( testClone );
  CPPUNIT_TEST_SUITE_END();

  // PUBLIC METHODS /////////////////////////////////////////////////////
public :

  /// Test clone()
  void testClone(void);

}; // class TestFaultCohesiveKin

#endif // pylith_faults_testfaultcohesivekin_hh


// End of file 
