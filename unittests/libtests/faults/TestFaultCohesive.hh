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
 * @file unittests/libtests/faults/TestFaultCohesive.hh
 *
 * @brief C++ TestFaultCohesive object
 *
 * C++ unit testing for Fault.
 */

#if !defined(pylith_faults_testfaultcohesive_hh)
#define pylith_faults_testfaultcohesive_hh

#include <cppunit/extensions/HelperMacros.h>

#include "pylith/utils/sievefwd.hh" // USES PETSc Mesh

/// Namespace for pylith package
namespace pylith {
  namespace faults {
    class Fault;
    class TestFaultCohesive;
    class CohesiveData;
  } // faults
} // pylith

/// C++ unit testing for Fault
class pylith::faults::TestFaultCohesive : public CppUnit::TestFixture
{ // class TestFaultCohesive

  // CPPUNIT TEST SUITE /////////////////////////////////////////////////
  CPPUNIT_TEST_SUITE( TestFaultCohesive );
  CPPUNIT_TEST( testAdjustTopologyLine2 );
  CPPUNIT_TEST( testAdjustTopologyTri3 );
  CPPUNIT_TEST( testAdjustTopologyTet4 );
  CPPUNIT_TEST_SUITE_END();

  // PUBLIC METHODS /////////////////////////////////////////////////////
public :

  /// Test adjustTopology() with 1-D line element.
  void testAdjustTopologyLine2(void);

  /// Test adjustTopology() with 2-D triangular element.
  void testAdjustTopologyTri3(void);

  /// Test adjustTopology() with 3-D tetrahedral element.
  void testAdjustTopologyTet4(void);

  // PROTECTED METHODS //////////////////////////////////////////////////
public :

  /// Test adjustTopology().
  void _testAdjustTopology(const CohesiveData& data);

}; // class TestFaultCohesive

#endif // pylith_faults_testfaultcohesive_hh

// End of file 
