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
// Copyright (c) 2010-2011 University of California, Davis
//
// See COPYING for license information.
//
// ----------------------------------------------------------------------
//

/**
 * @file unittests/libtests/meshio/TestOutputSolnPoints.hh
 *
 * @brief C++ TestOutputSolnPoints object
 *
 * C++ unit testing for OutputSolnPoints.
 */

#if !defined(pylith_meshio_testoutputsolnpoints_hh)
#define pylith_meshio_testoutputsolnpoints_hh

#include <cppunit/extensions/HelperMacros.h>

/// Namespace for pylith package
namespace pylith {
  namespace meshio {
    class TestOutputSolnPoints;
  } // meshio
} // pylith

/// C++ unit testing for OutputSolnPoints
class pylith::meshio::TestOutputSolnPoints : public CppUnit::TestFixture
{ // class TestOutputSolnPoints

  // CPPUNIT TEST SUITE /////////////////////////////////////////////////
  CPPUNIT_TEST_SUITE( TestOutputSolnPoints );

  CPPUNIT_TEST( testConstructor );
  CPPUNIT_TEST( testSetupInterpolator2D );
  CPPUNIT_TEST( testSetupInterpolator3D );

  CPPUNIT_TEST_SUITE_END();

  // PUBLIC METHODS /////////////////////////////////////////////////////
public :

  /// Test constructor
  void testConstructor(void);

  /// Test setupInterpolator for 2D mesh()
  void testSetupInterpolator2D(void);

  /// Test setupInterpolator for 3D mesh()
  void testSetupInterpolator3D(void);

}; // class TestOutputSolnPoints

#endif // pylith_meshio_testoutputsolnpoints_hh

// End of file 
