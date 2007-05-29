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

/// Namespace for pylith package
namespace pylith {
  namespace feassemble {
    class Integrator;
    class TestIntegrator;
    class IntegratorData;
  } // feassemble
} // pylith

/// C++ unit testing for Integrator
class pylith::feassemble::TestIntegrator : public CppUnit::TestFixture
{ // class TestIntegrator

  // CPPUNIT TEST SUITE /////////////////////////////////////////////////
  CPPUNIT_TEST_SUITE( TestIntegrator );

  CPPUNIT_TEST( testQuadrature );
  CPPUNIT_TEST( testNeedNewJacobian );
  CPPUNIT_TEST( testInitCellVector );
  CPPUNIT_TEST( testResetCellVector );
  CPPUNIT_TEST( testInitCellMatrix );
  CPPUNIT_TEST( testResetCellMatrix );

  CPPUNIT_TEST_SUITE_END();

  // PUBLIC METHODS /////////////////////////////////////////////////////
public :

  /// Test quadrature()
  void testQuadrature(void);

  /// Test needNewJacobian()
  void testNeedNewJacobian(void);

  /// Test _initCellVector()
  void testInitCellVector(void);

  /// Test _resetCellVector()
  void testResetCellVector(void);

  /// Test _initCellMatrix()
  void testInitCellMatrix(void);

  /// Test _resetCellMatrix()
  void testResetCellMatrix(void);

}; // class TestIntegrator

#endif // pylith_feassemble_testintegrator_hh

// End of file 
