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
    class FaultData;
  } // faults
} // pylith

/// C++ unit testing for Fault
class pylith::faults::TestFaultCohesive : public CppUnit::TestFixture
{ // class TestFaultCohesive

  // CPPUNIT TEST SUITE /////////////////////////////////////////////////
  CPPUNIT_TEST_SUITE( TestFaultCohesive );
  CPPUNIT_TEST( testAdjustTopologyLine );
  CPPUNIT_TEST_SUITE_END();

  // PUBLIC METHODS /////////////////////////////////////////////////////
public :

  /// Test adjustTopology() with 1-D line element.
  void testAdjustTopologyLine(void);

  // PROTECTED METHODS //////////////////////////////////////////////////
public :

  /// Test adjustTopology().
  void _testAdjustTopologyLine(const FaultData& data);

  /// Create mesh.
  void _createMesh(ALE::Obj<ALE::Mesh>* mesh,
		   const FaultData& data);

}; // class TestFaultCohesive

#endif // pylith_faults_testfaultcohesive_hh

// End of file 
