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

#include "TestIntegratorInertia.hh"

/// Namespace for spatialdata package
namespace pylith {
  namespace feassemble {
    class TestIntegratorInertia3D;
  } // feassemble
} // pylith

/// C++ unit testing for IntegratorInertia3D
class pylith::feassemble::TestIntegratorInertia3D : 
  public TestIntegratorInertia
{ // class TestIntegratorInertia3D

  // CPPUNIT TEST SUITE /////////////////////////////////////////////////
  CPPUNIT_TEST_SUITE( TestIntegratorInertia3D );

  CPPUNIT_TEST( testActionLinear );
  CPPUNIT_TEST( testIntegrateLinear );
  CPPUNIT_TEST( testLumpedLinear );

  CPPUNIT_TEST( testActionQuadratic );
  CPPUNIT_TEST( testIntegrateQuadratic );
  CPPUNIT_TEST( testLumpedQuadratic );

  CPPUNIT_TEST_SUITE_END();

  // PUBLIC METHODS /////////////////////////////////////////////////////
public :

  /// Test integrateAction() w/linear basis fns
  void testActionLinear(void);

  /// Test integrate() w/linear basis fns
  void testIntegrateLinear(void);

  /// Test integrateLumped() w/linear basis fns
  void testLumpedLinear(void);

  /// Test integrateAction() w/quadratic basis fns
  void testActionQuadratic(void);

  /// Test integrate() w/quadratic basis fns
  void testIntegrateQuadratic(void);

  /// Test integrateLumped() w/quadratic basis fns
  void testLumpedQuadratic(void);

}; // class TestIntegratorInertia3D

#endif // pylith_feassemble_testintegratorinertia3d_hh

// End of file 
