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
// Copyright (c) 2010-2012 University of California, Davis
//
// See COPYING for license information.
//
// ----------------------------------------------------------------------
//

/**
 * @file unittests/libtests/feassemble/TestQuadrature0D.hh
 *
 * @brief C++ TestQuadrature0D object
 *
 * C++ unit testing for Quadrature0D.
 */

#if !defined(pylith_feassemble_testquadrature0d_hh)
#define pylith_feassemble_testquadrature0d_hh

#include "TestQuadratureEngine.hh"

/// Namespace for pylith package
namespace pylith {
  namespace feassemble {
    class TestQuadrature0D;
  } // feassemble
} // pylith

/// C++ unit testing for Quadrature0D
class pylith::feassemble::TestQuadrature0D : public TestQuadratureEngine
{ // class TestQuadrature0D

  // CPPUNIT TEST SUITE /////////////////////////////////////////////////
  CPPUNIT_TEST_SUITE( TestQuadrature0D );

  CPPUNIT_TEST( testConstructor );
  CPPUNIT_TEST( testPoint );

  CPPUNIT_TEST_SUITE_END();

  // PUBLIC METHODS /////////////////////////////////////////////////////
public :

  /// Test constructor
  void testConstructor(void);

  /// Test initialize() & computeGeometry().
  void testPoint(void);

}; // class TestQuadrature0D

#endif // pylith_feassemble_testquadrature0d_hh

// End of file 
