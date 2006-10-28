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
 * @file unittests/libtests/feassemble/TestIntegratorInertia3D.hh
 *
 * @brief C++ TestIntegratorInertia3D object
 *
 * C++ unit testing for IntegratorInertia.
 */

#if !defined(pylith_feassemble_testintegratorinertia3d_hh)
#define pylith_feassemble_testintegratorinertia3d_hh

#include "TestIntegrator.hh"

/// Namespace for spatialdata package
namespace pylith {
  namespace feassemble {
    class TestIntegratorInertia3D;
  } // feassemble
} // pylith

/// C++ unit testing for IntegratorInertia3D
class pylith::feassemble::TestIntegratorInertia3D : public TestIntegrator
{ // class TestIntegratorInertia3D

  // CPPUNIT TEST SUITE /////////////////////////////////////////////////
  CPPUNIT_TEST_SUITE( TestIntegratorInertia3D );
  CPPUNIT_TEST( testLinear );
  CPPUNIT_TEST( testQuadratic );
  CPPUNIT_TEST_SUITE_END();

  // PUBLIC METHODS /////////////////////////////////////////////////////
public :

  /// Test integrate() & integrateAction() w/linear basis fns
  void testLinear(void);

  /// Test integrate() & integrateAction() w/quadratic basis fns
  void testQuadratic(void);

}; // class TestIntegratorInertia3D

#endif // pylith_feassemble_testintegratorinertia3d_hh

// End of file 
