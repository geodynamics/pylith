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

#include "TestIntegratorInertia.hh"

/// Namespace for pylith package
namespace pylith {
  namespace feassemble {
    class TestIntegratorInertia2D;
  } // feassemble
} // pylith

/// C++ unit testing for IntegratorInertia2D
class pylith::feassemble::TestIntegratorInertia2D : 
  public TestIntegratorInertia
{ // class TestIntegratorInertia2D

  // CPPUNIT TEST SUITE /////////////////////////////////////////////////
  CPPUNIT_TEST_SUITE( TestIntegratorInertia2D );

  CPPUNIT_TEST( testActionOne );
  CPPUNIT_TEST( testIntegrateOne );
  CPPUNIT_TEST( testLumpedOne );

  CPPUNIT_TEST( testActionOverlap1 );
  CPPUNIT_TEST( testIntegrateOverlap1 );
  CPPUNIT_TEST( testLumpedOverlap1 );

  CPPUNIT_TEST( testActionOverlap2 );
  CPPUNIT_TEST( testIntegrateOverlap2 );
  CPPUNIT_TEST( testLumpedOverlap2 );

  CPPUNIT_TEST_SUITE_END();

  // PUBLIC METHODS /////////////////////////////////////////////////////
public :

  /// Test integrateAction() w/linear basis fns (1 cell)
  void testActionOne(void);

  /// Test integrate() w/linear basis fns (1 cell)
  void testIntegrateOne(void);

  /// Test integrateLumped() w/linear basis fns (1 cell)
  void testLumpedOne(void);

  /// Test integrateAction() w/linear basis fns (2 cells sharing 1 vertex)
  void testActionOverlap1(void);

  /// Test integrate() w/linear basis fns (2 cells sharing 1 vertex)
  void testIntegrateOverlap1(void);

  /// Test integrateLumped() w/linear basis fns (2 cells sharing 1 vertex)
  void testLumpedOverlap1(void);

  /// Test integrateAction() w/linear basis fns (2 cells sharing 2 vertices)
  void testActionOverlap2(void);

  /// Test integrate() w/linear basis fns (2 cells sharing 2 vertices)
  void testIntegrateOverlap2(void);

  /// Test integrateLumped() w/linear basis fns (2 cells sharing 2 vertices)
  void testLumpedOverlap2(void);

}; // class TestIntegratorInertia2D

#endif // pylith_feassemble_testintegratorinertia2d_hh

// End of file 
