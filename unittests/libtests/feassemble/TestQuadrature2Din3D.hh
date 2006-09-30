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
 * @file unittests/libtests/feassemble/TestQuadrature2Din3D.hh
 *
 * @brief C++ TestQuadrature2Din3D object
 *
 * C++ unit testing for Quadrature2Din3D.
 */

#if !defined(pylith_feassemble_testquadrature2din3d_hh)
#define pylith_feassemble_testquadrature2din3d_hh

#include "TestQuadrature.hh"

/// Namespace for spatialdata package
namespace pylith {
  namespace feassemble {
    class TestQuadrature2Din3D;
  } // feassemble
} // pylith

/// C++ unit testing for Quadrature2Din3D
class pylith::feassemble::TestQuadrature2Din3D : public TestQuadrature
{ // class TestQuadrature2D

  // CPPUNIT TEST SUITE /////////////////////////////////////////////////
  CPPUNIT_TEST_SUITE( TestQuadrature2Din3D );
  CPPUNIT_TEST( testConstructor );
  CPPUNIT_TEST( testLinearXYZ );
  CPPUNIT_TEST( testLinearXY );
  CPPUNIT_TEST( testLinearYZ );
  CPPUNIT_TEST( testLinearXZ );
  CPPUNIT_TEST( testQuadratic );
  CPPUNIT_TEST_SUITE_END();

  // PUBLIC METHODS /////////////////////////////////////////////////////
public :

  /// Test constructor
  void testConstructor(void);

  /// Test initialize() & computeGeometry() w/linear basis fns
  void testLinearXYZ(void);

  /// Test initialize() & computeGeometry() w/linear basis fns
  void testLinearXY(void);

  /// Test initialize() & computeGeometry() w/linear basis fns
  void testLinearYZ(void);

  /// Test initialize() & computeGeometry() w/linear basis fns
  void testLinearXZ(void);

  /// Test initialize() & computeGeometry() w/quadratic basis fns
  void testQuadratic(void);

}; // class TestQuadrature

#endif // pylith_feassemble_testquadrature2din3d_hh

// End of file 
