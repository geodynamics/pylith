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
 * @file unittests/libtests/faults/TestFaultCohesiveKinTet4e.hh
 *
 * @brief C++ TestFaultCohesiveKinTet4e object.
 *
 * C++ unit testing for FaultCohesiveKin for mesh with 3-D tetrahedral cells.
 */

#if !defined(pylith_faults_testfaultcohesivekintet4e_hh)
#define pylith_faults_testfaultcohesivekintet4e_hh

#include "TestFaultCohesiveKin.hh" // ISA TestFaultCohesiveKin

/// Namespace for pylith package
namespace pylith {
  namespace faults {
    class TestFaultCohesiveKinTet4e;
  } // bc
} // pylith

/// C++ unit testing for FaultCohesiveKin for mesh with 3-D tetrahedral cells.
class pylith::faults::TestFaultCohesiveKinTet4e : public TestFaultCohesiveKin
{ // class TestFaultCohesiveKinTet4e

  // CPPUNIT TEST SUITE /////////////////////////////////////////////////
  CPPUNIT_TEST_SUITE( TestFaultCohesiveKinTet4e );

  CPPUNIT_TEST( testInitialize );
  CPPUNIT_TEST( testIntegrateResidual );
  CPPUNIT_TEST( testIntegrateJacobian );

  CPPUNIT_TEST_SUITE_END();

  // PUBLIC METHODS /////////////////////////////////////////////////////
public :

  /// Setup testing data.
  void setUp(void);

}; // class TestFaultCohesiveKinTet4e

#endif // pylith_faults_testfaultcohesivetet4e_hh


// End of file 
