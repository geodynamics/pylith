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
 * @file unittests/libtests/feassemble/TestIntegratorInertia1D.hh
 *
 * @brief C++ TestIntegratorInertia1D object
 *
 * C++ unit testing for Quadrature1D.
 */

#if !defined(pylith_feassemble_testintegratorinertia1d_hh)
#define pylith_feassemble_testintegratorinertia1d_hh

#include "TestIntegrator.hh"

/// Namespace for spatialdata package
namespace pylith {
  namespace feassemble {
    class TestIntegratorInertia1D;
  } // feassemble
} // pylith

/// C++ unit testing for IntegratorInertia1D
class pylith::feassemble::TestIntegratorInertia1D : public TestIntegrator
{ // class TestIntegratorInertia1D

  // CPPUNIT TEST SUITE /////////////////////////////////////////////////
  CPPUNIT_TEST_SUITE( TestIntegratorInertia1D );
  CPPUNIT_TEST( testConstructor );
  CPPUNIT_TEST( testLinear );
  CPPUNIT_TEST( testQuadratic );
  CPPUNIT_TEST_SUITE_END();

  // PUBLIC METHODS /////////////////////////////////////////////////////
public :

  /// Test constructor
  void testConstructor(void);

  /// Test integrate() & integrateAction() w/linear basis fns
  void testLinear(void);

  /// Test integrate() & integrateAction() w/quadratic basis fns
  void testQuadratic(void);

}; // class TestIntegratorInertia1D

#endif // pylith_feassemble_testintegratorinertia1d_hh

// End of file 
