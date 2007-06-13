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
 * @file unittests/libtests/faults/TestFaultCohesiveKinHex8.hh
 *
 * @brief C++ TestFaultCohesiveKinHex8 object.
 *
 * C++ unit testing for FaultCohesiveKin for mesh with 3-D hex cells.
 */

#if !defined(pylith_faults_testfaultcohesivekinhex8_hh)
#define pylith_faults_testfaultcohesivekinhex8_hh

#include "TestFaultCohesiveKin.hh" // ISA TestFaultCohesiveKin

/// Namespace for pylith package
namespace pylith {
  namespace faults {
    class TestFaultCohesiveKinHex8;
  } // bc
} // pylith

/// C++ unit testing for FaultCohesiveKin for mesh with 3-D hex cells.
class pylith::faults::TestFaultCohesiveKinHex8 : public TestFaultCohesiveKin
{ // class TestFaultCohesiveKinHex8

  // CPPUNIT TEST SUITE /////////////////////////////////////////////////
  CPPUNIT_TEST_SUITE( TestFaultCohesiveKinHex8 );

  CPPUNIT_TEST( testInitialize );
  CPPUNIT_TEST( testIntegrateResidual );
  CPPUNIT_TEST( testIntegrateJacobian );

  CPPUNIT_TEST_SUITE_END();

  // PUBLIC METHODS /////////////////////////////////////////////////////
public :

  /// Setup testing data.
  void setUp(void);

}; // class TestFaultCohesiveKinHex8

#endif // pylith_faults_testfaultcohesivehex8_hh


// End of file 
