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
 * @file unittests/libtests/faults/TestFaultCohesiveKinTet4f.hh
 *
 * @brief C++ TestFaultCohesiveKinTet4f object.
 *
 * C++ unit testing for FaultCohesiveKin for mesh with 3-D tetrahedral cells.
 */

#if !defined(pylith_faults_testfaultcohesivekintet4f_hh)
#define pylith_faults_testfaultcohesivekintet4f_hh

#include "TestFaultCohesiveKin.hh" // ISA TestFaultCohesiveKin

/// Namespace for pylith package
namespace pylith {
  namespace faults {
    class TestFaultCohesiveKinTet4f;
  } // bc
} // pylith

/// C++ unit testing for FaultCohesiveKin for mesh with 3-D tetrahedral cells.
class pylith::faults::TestFaultCohesiveKinTet4f : public TestFaultCohesiveKin
{ // class TestFaultCohesiveKinTet4f

  // CPPUNIT TEST SUITE /////////////////////////////////////////////////
  CPPUNIT_TEST_SUITE( TestFaultCohesiveKinTet4f );

  CPPUNIT_TEST( testInitialize );
  CPPUNIT_TEST( testIntegrateResidual );
  CPPUNIT_TEST( testIntegrateJacobian );

  CPPUNIT_TEST_SUITE_END();

  // PUBLIC METHODS /////////////////////////////////////////////////////
public :

  /// Setup testing data.
  void setUp(void);

}; // class TestFaultCohesiveKinTet4f

#endif // pylith_faults_testfaultcohesivetet4f_hh


// End of file 
