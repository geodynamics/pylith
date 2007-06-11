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
 * @file unittests/libtests/faults/TestFaultCohesiveKinTet4.hh
 *
 * @brief C++ TestFaultCohesiveKinTet4 object.
 *
 * C++ unit testing for FaultCohesiveKin for mesh with 3-D tetrahedral cells.
 */

#if !defined(pylith_faults_testfaultcohesivekintet4_hh)
#define pylith_faults_testfaultcohesivekintet4_hh

#include "TestFaultCohesiveKin.hh" // ISA TestFaultCohesiveKin

/// Namespace for pylith package
namespace pylith {
  namespace faults {
    class TestFaultCohesiveKinTet4;
  } // bc
} // pylith

/// C++ unit testing for FaultCohesiveKin for mesh with 3-D tetrahedral cells.
class pylith::faults::TestFaultCohesiveKinTet4 : public TestFaultCohesiveKin
{ // class TestFaultCohesiveKinTet4

  // CPPUNIT TEST SUITE /////////////////////////////////////////////////
  CPPUNIT_TEST_SUITE( TestFaultCohesiveKinTet4 );

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

}; // class TestFaultCohesiveKinTet4

#endif // pylith_faults_testfaultcohesivetet4_hh


// End of file 
