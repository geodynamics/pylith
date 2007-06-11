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
 * @file unittests/libtests/faults/TestFaultCohesiveKinQuad4.hh
 *
 * @brief C++ TestFaultCohesiveKinQuad4 object.
 *
 * C++ unit testing for FaultCohesiveKin for mesh with 2-D quadrilateral cells.
 */

#if !defined(pylith_faults_testfaultcohesivekinquad4_hh)
#define pylith_faults_testfaultcohesivekinquad4_hh

#include "TestFaultCohesiveKin.hh" // ISA TestFaultCohesiveKin

/// Namespace for pylith package
namespace pylith {
  namespace faults {
    class TestFaultCohesiveKinQuad4;
  } // bc
} // pylith

/// C++ unit testing for FaultCohesiveKin for mesh with 2-D quadrilateral cells.
class pylith::faults::TestFaultCohesiveKinQuad4 : public TestFaultCohesiveKin
{ // class TestFaultCohesiveKinQuad4

  // CPPUNIT TEST SUITE /////////////////////////////////////////////////
  CPPUNIT_TEST_SUITE( TestFaultCohesiveKinQuad4 );

  CPPUNIT_TEST( testInitialize );
  CPPUNIT_TEST( testIntegrateResidual );
  CPPUNIT_TEST( testIntegrateJacobian );
  CPPUNIT_TEST( testSetConstraintSizes );
  CPPUNIT_TEST( testSetField );

  CPPUNIT_TEST_SUITE_END();

  // PUBLIC METHODS /////////////////////////////////////////////////////
public :

  /// Setup testing data.
  void setUp(void);

}; // class TestFaultCohesiveKinQuad4

#endif // pylith_faults_testfaultcohesivequad4_hh


// End of file 
