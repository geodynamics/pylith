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
  CPPUNIT_TEST( testAdjustTopologyTri3b );
  CPPUNIT_TEST( testAdjustTopologyTri3c );
  CPPUNIT_TEST( testAdjustTopologyQuad4 );
  CPPUNIT_TEST( testAdjustTopologyTet4 );
  CPPUNIT_TEST( testAdjustTopologyHex8 );

  CPPUNIT_TEST( testAdjustTopologyLine2Lagrange );
  CPPUNIT_TEST( testAdjustTopologyTri3Lagrange );
  CPPUNIT_TEST( testAdjustTopologyQuad4Lagrange );
  CPPUNIT_TEST( testAdjustTopologyTet4Lagrange );
  CPPUNIT_TEST( testAdjustTopologyHex8Lagrange );

  CPPUNIT_TEST( testOrient1D );
  CPPUNIT_TEST( testOrient2D );
  CPPUNIT_TEST( testOrient3D );

  CPPUNIT_TEST_SUITE_END();

  // PUBLIC METHODS /////////////////////////////////////////////////////
public :

  /// Test adjustTopology() with 1-D line element.
  void testAdjustTopologyLine2(void);

  /// Test adjustTopology() with 2-D triangular element.
  void testAdjustTopologyTri3(void);

  /// Test adjustTopology() with 2-D triangular element (edge b).
  void testAdjustTopologyTri3b(void);

  /// Test adjustTopology() with 2-D triangular element (edge c).
  void testAdjustTopologyTri3c(void);

  /// Test adjustTopology() with 2-D quadrilateral element.
  void testAdjustTopologyQuad4(void);

  /// Test adjustTopology() with 3-D tetrahedral element.
  void testAdjustTopologyTet4(void);

  /// Test adjustTopology() with 3-D hexahedral element.
  void testAdjustTopologyHex8(void);

  /// Test adjustTopology() with 1-D line element for Lagrange
  /// multipliers.
  void testAdjustTopologyLine2Lagrange(void);

  /// Test adjustTopology() with 2-D triangular element for Lagrange
  /// multipliers.
  void testAdjustTopologyTri3Lagrange(void);

  /// Test adjustTopology() with 2-D quadrilateral element for Lagrange
  /// multipliers.
  void testAdjustTopologyQuad4Lagrange(void);

  /// Test adjustTopology() with 3-D tetrahedral element for Lagrange
  /// multipliers.
  void testAdjustTopologyTet4Lagrange(void);

  /// Test adjustTopology() with 3-D hexahedral element for Lagrange
  /// multipliers.
  void testAdjustTopologyHex8Lagrange(void);

  /// Test _orientationSize().
  void testOrientationSize(void);

  /// Test _orient1D().
  void testOrient1D(void);

  /// Test _orient2D().
  void testOrient2D(void);

  /// Test _orient3D().
  void testOrient3D(void);

  // PROTECTED METHODS //////////////////////////////////////////////////
public :

  /** Test adjustTopology().
   *
   * @param fault Fault for cohesive elements.
   * @param data Cohesive element data.
   */
  void _testAdjustTopology(Fault* fault,
			   const CohesiveData& data);

}; // class TestFaultCohesive

#endif // pylith_faults_testfaultcohesive_hh

// End of file 
