// -*- C++ -*-
//
// ----------------------------------------------------------------------
//
// Brad T. Aagaard, U.S. Geological Survey
// Charles A. Williams, GNS Science
// Matthew G. Knepley, University at Buffalo
//
// This code was developed as part of the Computational Infrastructure
// for Geodynamics (http://geodynamics.org).
//
// Copyright (c) 2010-2021 University of California, Davis
//
// See LICENSE.md for license information.
//
// ----------------------------------------------------------------------
//

/**
 * @file tests/libtests/meshio/TestVertexFilterVecNorm.hh
 *
 * @brief C++ TestVertexFilterVecNorm object
 *
 * C++ unit testing for VertexFilterVecNorm.
 */

#if !defined(pylith_meshio_testvertexfiltervecnorm_hh)
#define pylith_meshio_testvertexfiltervecnorm_hh

#include <cppunit/extensions/HelperMacros.h>

/// Namespace for pylith package
namespace pylith {
  namespace meshio {
    class TestVertexFilterVecNorm;
  } // meshio
} // pylith

/// C++ unit testing for VertexFilterVecNorm
class pylith::meshio::TestVertexFilterVecNorm : public CppUnit::TestFixture
{ // class TestVertexFilterVecNorm

  // CPPUNIT TEST SUITE /////////////////////////////////////////////////
  CPPUNIT_TEST_SUITE( TestVertexFilterVecNorm );

  CPPUNIT_TEST( testConstructor );
  CPPUNIT_TEST( testFilter );

  CPPUNIT_TEST_SUITE_END();

  // PUBLIC METHODS /////////////////////////////////////////////////////
public :

  /// Test constructor
  void testConstructor(void);

  /// Test filter()
  void testFilter(void);

}; // class TestVertexFilterVecNorm

#endif // pylith_meshio_testvertexfiltervecnorm_hh

// End of file 
