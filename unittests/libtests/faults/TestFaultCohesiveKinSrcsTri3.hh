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
 * @file unittests/libtests/faults/TestFaultCohesiveKinSrcsTri3.hh
 *
 * @brief C++ TestFaultCohesiveKinSrcsTri3 object.
 *
 * C++ unit testing for FaultCohesiveKinSrcs for mesh with 2-D triangular cells.
 */

#if !defined(pylith_faults_testfaultcohesivekinsrcstri3_hh)
#define pylith_faults_testfaultcohesivekinsrcstri3_hh

#include "TestFaultCohesiveKinSrcs.hh" // ISA TestFaultCohesiveKinSrcs

/// Namespace for pylith package
namespace pylith {
  namespace faults {
    class TestFaultCohesiveKinSrcsTri3;
  } // bc
} // pylith

/// C++ unit testing for FaultCohesiveKinSrcs for mesh with 2-D triangular cells.
class pylith::faults::TestFaultCohesiveKinSrcsTri3 : public TestFaultCohesiveKinSrcs
{ // class TestFaultCohesiveKinSrcsTri3

  // CPPUNIT TEST SUITE /////////////////////////////////////////////////
  CPPUNIT_TEST_SUITE( TestFaultCohesiveKinSrcsTri3 );

  CPPUNIT_TEST( testInitialize );
  CPPUNIT_TEST( testIntegrateResidual );
  CPPUNIT_TEST( testIntegrateJacobian );
  CPPUNIT_TEST( testCalcTractionsChange );

  CPPUNIT_TEST_SUITE_END();

  // PUBLIC METHODS /////////////////////////////////////////////////////
public :

  /// Setup testing data.
  void setUp(void);

}; // class TestFaultCohesiveKinSrcsTri3

#endif // pylith_faults_testfaultcohesivekinsrcstri3_hh


// End of file 
