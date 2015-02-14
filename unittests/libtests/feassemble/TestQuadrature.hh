// -*- C++ -*-
//
// ----------------------------------------------------------------------
//
// Brad T. Aagaard, U.S. Geological Survey
// Charles A. Williams, GNS Science
// Matthew G. Knepley, University of Chicago
//
// This code was developed as part of the Computational Infrastructure
// for Geodynamics (http://geodynamics.org).
//
// Copyright (c) 2010-2015 University of California, Davis
//
// See COPYING for license information.
//
// ----------------------------------------------------------------------
//

/**
 * @file unittests/libtests/feassemble/TestQuadrature.hh
 *
 * @brief C++ TestQuadrature object
 *
 * C++ unit testing for Quadrature.
 */

#if !defined(pylith_feassemble_testquadrature_hh)
#define pylith_feassemble_testquadrature_hh

#include <cppunit/extensions/HelperMacros.h>

/// Namespace for pylith package
namespace pylith {
  namespace feassemble {
    class TestQuadrature;
  } // feassemble
} // pylith

/// C++ unit testing for Quadrature
class pylith::feassemble::TestQuadrature : public CppUnit::TestFixture
{ // class TestQuadrature

  // CPPUNIT TEST SUITE /////////////////////////////////////////////////
  CPPUNIT_TEST_SUITE( TestQuadrature );

  CPPUNIT_TEST( testCopyConstructor );
  CPPUNIT_TEST( testCheckConditioning );
  CPPUNIT_TEST( testEngineAccessors );
  CPPUNIT_TEST( testComputeGeometryCell );

  CPPUNIT_TEST_SUITE_END();

  // PUBLIC METHODS /////////////////////////////////////////////////////
public :

  /// Test copy constructor.
  void testCopyConstructor(void);

  /// Test checkConditioning()
  void testCheckConditioning(void);

  /// Test quadPts(), basisDeriv(), jacobian(), and jacobianDet().
  void testEngineAccessors(void);

  /// Test computeGeometry() with coordinates and cell.
  void testComputeGeometryCell(void);

}; // class TestQuadrature

#endif // pylith_feassemble_testquadrature_hh

// End of file 
