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
 * @file unittests/libtests/faults/TestFaultCohesiveKinLine2.hh
 *
 * @brief C++ TestFaultCohesiveKinLine2 object.
 *
 * C++ unit testing for FaultCohesiveKin for mesh with 1-D line cells.
 */

#if !defined(pylith_faults_testfaultcohesivekinline2_hh)
#define pylith_faults_testfaultcohesivekinline2_hh

#include "TestFaultCohesiveKin.hh" // ISA TestFaultCohesiveKin

/// Namespace for pylith package
namespace pylith {
  namespace faults {
    class TestFaultCohesiveKinLine2;
  } // bc
} // pylith

/// C++ unit testing for FaultCohesiveKin for mesh with 1-D line cells.
class pylith::faults::TestFaultCohesiveKinLine2 : public TestFaultCohesiveKin
{ // class TestFaultCohesiveKinLine2

  // CPPUNIT TEST SUITE /////////////////////////////////////////////////
  CPPUNIT_TEST_SUITE( TestFaultCohesiveKinLine2 );

  CPPUNIT_TEST( testInitialize );
  CPPUNIT_TEST( testIntegrateResidual );
  CPPUNIT_TEST( testIntegrateJacobian );
  CPPUNIT_TEST( testSetConstraintSizes );
  CPPUNIT_TEST( testSetConstraints );
  CPPUNIT_TEST( testSetField );

  CPPUNIT_TEST_SUITE_END();

  // PUBLIC METHODS /////////////////////////////////////////////////////
public :

  /// Setup testing data.
  void setUp(void);

  /// Tear down testing data.
  void tearDown(void);

}; // class TestFaultCohesiveKinLine2

#endif // pylith_faults_testfaultcohesiveline2_hh


// End of file 
