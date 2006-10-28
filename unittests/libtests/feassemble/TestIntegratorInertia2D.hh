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
 * @file unittests/libtests/feassemble/TestIntegratorInertia2D.hh
 *
 * @brief C++ TestIntegratorInertia2D object
 *
 * C++ unit testing for IntegratorInertia.
 */

#if !defined(pylith_feassemble_testintegratorinertia2d_hh)
#define pylith_feassemble_testintegratorinertia2d_hh

#include "TestIntegrator.hh"

/// Namespace for spatialdata package
namespace pylith {
  namespace feassemble {
    class TestIntegratorInertia2D;
  } // feassemble
} // pylith

/// C++ unit testing for IntegratorInertia2D
class pylith::feassemble::TestIntegratorInertia2D : public TestIntegrator
{ // class TestIntegratorInertia2D

  // CPPUNIT TEST SUITE /////////////////////////////////////////////////
  CPPUNIT_TEST_SUITE( TestIntegratorInertia2D );
  CPPUNIT_TEST( testOne );
  CPPUNIT_TEST( testOverlap1 );
  CPPUNIT_TEST( testOverlap2 );
  CPPUNIT_TEST_SUITE_END();

  // PUBLIC METHODS /////////////////////////////////////////////////////
public :

  /// Test integrate() & integrateAction() w/linear basis fns
  /// (1 cell)
  void testOne(void);

  /// Test integrate() & integrateAction() w/linear basis fns
  /// (2 cells sharing 1 vertex)
  void testOverlap1(void);

  /// Test integrate() & integrateAction() w/linear basis fns
  /// (2 cells sharing 2 vertices)
  void testOverlap2(void);

}; // class TestIntegratorInertia2D

#endif // pylith_feassemble_testintegratorinertia2d_hh

// End of file 
