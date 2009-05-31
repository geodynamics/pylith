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
 * @file unittests/libtests/feassemble/TestIntegrator.hh
 *
 * @brief C++ TestIntegrator object
 *
 * C++ unit testing for Integrator.
 */

#if !defined(pylith_feassemble_testintegrator_hh)
#define pylith_feassemble_testintegrator_hh

#include <cppunit/extensions/HelperMacros.h>

#include "pylith/topology/topologyfwd.hh" // USES Mesh
#include "pylith/feassemble/feassemblefwd.hh" // USES Quadrature

/// Namespace for pylith package
namespace pylith {
  namespace feassemble {
    class TestIntegrator;
  } // feassemble
} // pylith

/// C++ unit testing for Integrator
class pylith::feassemble::TestIntegrator : public CppUnit::TestFixture
{ // class TestIntegrator

  // CPPUNIT TEST SUITE /////////////////////////////////////////////////
  CPPUNIT_TEST_SUITE( TestIntegrator );

  CPPUNIT_TEST( testTimeStep );
  CPPUNIT_TEST( testStableTimeStep );  

  CPPUNIT_TEST( testQuadrature );
  CPPUNIT_TEST( testNormalizer );
  CPPUNIT_TEST( testGravityField );
  CPPUNIT_TEST( testInitCellVector );
  CPPUNIT_TEST( testResetCellVector );
  CPPUNIT_TEST( testInitCellMatrix );
  CPPUNIT_TEST( testResetCellMatrix );
  CPPUNIT_TEST( testSplitFields );

  CPPUNIT_TEST_SUITE_END();

  // PUBLIC METHODS /////////////////////////////////////////////////////
public :

  /// Test timeStep().
  void testTimeStep(void);

  /// Test stableTimeStep().
  void testStableTimeStep(void);

  /// Test quadrature().
  void testQuadrature(void);

  /// Test normalizer().
  void testNormalizer(void);

  /// Test gravityField().
  void testGravityField(void);

  /// Test _initCellVector().
  void testInitCellVector(void);

  /// Test _resetCellVector().
  void testResetCellVector(void);

  /// Test _initCellMatrix().
  void testInitCellMatrix(void);

  /// Test _resetCellMatrix().
  void testResetCellMatrix(void);

  /// Test splitFields().
  void testSplitFields(void);

  // PRIVATE METHODS /////////////////////////////////////////////////////
private :

  /** Initialize quadrature.
   *
   * @param quadrature Quadrature to initiqlize.
   */
  void _initQuadrature(Quadrature<topology::Mesh>* quadrature);

}; // class TestIntegrator

#endif // pylith_feassemble_testintegrator_hh

// End of file 
