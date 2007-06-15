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
 * @file unittests/libtests/faults/TestFaultCohesiveKinQuad4e.hh
 *
 * @brief C++ TestFaultCohesiveKinQuad4e object.
 *
 * C++ unit testing for FaultCohesiveKin for mesh with 2-D quadrilateral cells.
 */

#if !defined(pylith_faults_testfaultcohesivekinquad4e_hh)
#define pylith_faults_testfaultcohesivekinquad4e_hh

#include "TestFaultCohesiveKin.hh" // ISA TestFaultCohesiveKin

/// Namespace for pylith package
namespace pylith {
  namespace faults {
    class TestFaultCohesiveKinQuad4e;
  } // bc
} // pylith

/// C++ unit testing for FaultCohesiveKin for mesh with 2-D quadrilateral cells.
class pylith::faults::TestFaultCohesiveKinQuad4e : public TestFaultCohesiveKin
{ // class TestFaultCohesiveKinQuad4e

  // CPPUNIT TEST SUITE /////////////////////////////////////////////////
  CPPUNIT_TEST_SUITE( TestFaultCohesiveKinQuad4e );

  CPPUNIT_TEST( testInitialize );
  CPPUNIT_TEST( testIntegrateResidual );
  CPPUNIT_TEST( testIntegrateJacobian );

  CPPUNIT_TEST_SUITE_END();

  // PUBLIC METHODS /////////////////////////////////////////////////////
public :

  /// Setup testing data.
  void setUp(void);

}; // class TestFaultCohesiveKinQuad4e

#endif // pylith_faults_testfaultcohesivequad4e_hh


// End of file 
