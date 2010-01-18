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
 * @file unittests/libtests/faults/TestFaultCohesiveDynLQuad4.hh
 *
 * @brief C++ TestFaultCohesiveDynLQuad4 object.
 *
 * C++ unit testing for FaultCohesiveDynL for mesh with 2-D quadrilateral cells.
 */

#if !defined(pylith_faults_testfaultcohesivedynlquad4_hh)
#define pylith_faults_testfaultcohesivedynlquad4_hh

#include "TestFaultCohesiveDynL.hh" // ISA TestFaultCohesiveDynL

/// Namespace for pylith package
namespace pylith {
  namespace faults {
    class TestFaultCohesiveDynLQuad4;
  } // bc
} // pylith

/// C++ unit testing for FaultCohesiveDynL for mesh with 2-D quadrilateral cells.
class pylith::faults::TestFaultCohesiveDynLQuad4 : public TestFaultCohesiveDynL
{ // class TestFaultCohesiveDynLQuad4

  // CPPUNIT TEST SUITE /////////////////////////////////////////////////
  CPPUNIT_TEST_SUITE( TestFaultCohesiveDynLQuad4 );

  CPPUNIT_TEST( testInitialize );
  CPPUNIT_TEST( testConstrainSolnSpaceStick );
  CPPUNIT_TEST( testConstrainSolnSpaceSlip );
  CPPUNIT_TEST( testConstrainSolnSpaceOpen );
  CPPUNIT_TEST( testUpdateStateVars );
  CPPUNIT_TEST( testCalcTractions );

  CPPUNIT_TEST_SUITE_END();

  // PUBLIC METHODS /////////////////////////////////////////////////////
public :

  /// Setup testing data.
  void setUp(void);

}; // class TestFaultCohesiveDynLQuad4

#endif // pylith_faults_testfaultcohesivedynlquad4_hh


// End of file 
