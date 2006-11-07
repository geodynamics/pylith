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
 * C++ unit testing for IntegratorInertia.
 */

#if !defined(pylith_feassemble_testintegratorinertia1d_hh)
#define pylith_feassemble_testintegratorinertia1d_hh

#include "TestIntegratorInertia.hh"

/// Namespace for pylith package
namespace pylith {
  namespace feassemble {
    class TestIntegratorInertia1D;
  } // feassemble
} // pylith

/// C++ unit testing for IntegratorInertia1D
class pylith::feassemble::TestIntegratorInertia1D : 
  public TestIntegratorInertia
{ // class TestIntegratorInertia1D

  // CPPUNIT TEST SUITE /////////////////////////////////////////////////
  CPPUNIT_TEST_SUITE( TestIntegratorInertia1D );
  CPPUNIT_TEST( testConstructor );
  CPPUNIT_TEST( testActionLinear );
  CPPUNIT_TEST( testIntegrateLinear );
  CPPUNIT_TEST( testLumpedLinear );
  CPPUNIT_TEST( testActionQuadratic );
  CPPUNIT_TEST( testIntegrateQuadratic );
  CPPUNIT_TEST( testLumpedQuadratic );
  CPPUNIT_TEST_SUITE_END();

  // PUBLIC METHODS /////////////////////////////////////////////////////
public :

  /// Test constructor
  void testConstructor(void);

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

}; // class TestIntegratorInertia1D

#endif // pylith_feassemble_testintegratorinertia1d_hh

// End of file 
