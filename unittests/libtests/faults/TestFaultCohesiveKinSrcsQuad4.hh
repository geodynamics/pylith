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
 * @file unittests/libtests/faults/TestFaultCohesiveKinSrcsQuad4.hh
 *
 * @brief C++ TestFaultCohesiveKinSrcsQuad4 object.
 *
 * C++ unit testing for FaultCohesiveKin for mesh with 2-D quadrilateral cells.
 */

#if !defined(pylith_faults_testfaultcohesivekinsrcsquad4_hh)
#define pylith_faults_testfaultcohesivekinsrcsquad4_hh

#include "TestFaultCohesiveKinSrcs.hh" // ISA TestFaultCohesiveKin

/// Namespace for pylith package
namespace pylith {
  namespace faults {
    class TestFaultCohesiveKinSrcsQuad4;
  } // bc
} // pylith

/// C++ unit testing for FaultCohesiveKin for mesh with 2-D quadrilateral cells.
class pylith::faults::TestFaultCohesiveKinSrcsQuad4 : public TestFaultCohesiveKinSrcs
{ // class TestFaultCohesiveKinSrcsQuad4

  // CPPUNIT TEST SUITE /////////////////////////////////////////////////
  CPPUNIT_TEST_SUITE( TestFaultCohesiveKinSrcsQuad4 );

  CPPUNIT_TEST( testInitialize );
  CPPUNIT_TEST( testIntegrateResidual );
  CPPUNIT_TEST( testIntegrateJacobian );
  CPPUNIT_TEST( testCalcTractionsChange );

  CPPUNIT_TEST_SUITE_END();

  // PUBLIC METHODS /////////////////////////////////////////////////////
public :

  /// Setup testing data.
  void setUp(void);

}; // class TestFaultCohesiveKinSrcsQuad4

#endif // pylith_faults_testfaultcohesivekinsrcsquad4_hh


// End of file 
