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
 * @file unittests/libtests/feassemble/TestQuadrature1Din3D.hh
 *
 * @brief C++ TestQuadrature1Din3D object
 *
 * C++ unit testing for Quadrature1Din3D.
 */

#if !defined(pylith_feassemble_testquadrature1din3d_hh)
#define pylith_feassemble_testquadrature1din3d_hh

#include "TestQuadrature.hh"

/// Namespace for pylith package
namespace pylith {
  namespace feassemble {
    class TestQuadrature1Din3D;
  } // feassemble
} // pylith

/// C++ unit testing for Quadrature1Din3D
class pylith::feassemble::TestQuadrature1Din3D : public TestQuadrature
{ // class TestQuadrature1D

  // CPPUNIT TEST SUITE /////////////////////////////////////////////////
  CPPUNIT_TEST_SUITE( TestQuadrature1Din3D );
  CPPUNIT_TEST( testConstructor );
  CPPUNIT_TEST( testLinear );
  CPPUNIT_TEST( testQuadratic );
  CPPUNIT_TEST_SUITE_END();

  // PUBLIC METHODS /////////////////////////////////////////////////////
public :

  /// Test constructor
  void testConstructor(void);

  /// Test initialize() & computeGeometry() w/linear basis fns
  void testLinear(void);

  /// Test initialize() & computeGeometry() w/quadratic basis fns
  void testQuadratic(void);

}; // class TestQuadrature

#endif // pylith_feassemble_testquadrature1din3d_hh

// End of file 
