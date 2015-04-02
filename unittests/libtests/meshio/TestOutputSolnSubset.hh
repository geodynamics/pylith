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
 * @file unittests/libtests/meshio/TestOutputSolnSubset.hh
 *
 * @brief C++ TestOutputSolnSubset object
 *
 * C++ unit testing for OutputSolnSubset.
 */

#if !defined(pylith_meshio_testoutputsolnsubset_hh)
#define pylith_meshio_testoutputsolnsubset_hh

#include <cppunit/extensions/HelperMacros.h>

/// Namespace for pylith package
namespace pylith {
  namespace meshio {
    class TestOutputSolnSubset;
  } // meshio
} // pylith

/// C++ unit testing for OutputSolnSubset
class pylith::meshio::TestOutputSolnSubset : public CppUnit::TestFixture
{ // class TestOutputSolnSubset

  // CPPUNIT TEST SUITE /////////////////////////////////////////////////
  CPPUNIT_TEST_SUITE( TestOutputSolnSubset );

  CPPUNIT_TEST( testConstructor );
  CPPUNIT_TEST( testLabel );
  CPPUNIT_TEST( testSubdomainMesh );

  CPPUNIT_TEST_SUITE_END();

  // PUBLIC METHODS /////////////////////////////////////////////////////
public :

  /// Test constructor
  void testConstructor(void);

  /// Test label()
  void testLabel(void);

  /// Test subdomainMesh()
  void testSubdomainMesh(void);

}; // class TestOutputSolnSubset

#endif // pylith_meshio_testoutputsolnsubset_hh

// End of file 
