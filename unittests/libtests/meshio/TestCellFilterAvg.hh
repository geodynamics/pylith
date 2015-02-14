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
 * @file unittests/libtests/meshio/TestCellFilterAvg.hh
 *
 * @brief C++ TestCellFilterAvg object
 *
 * C++ unit testing for CellFilterAvg.
 */

#if !defined(pylith_meshio_testcellfilteravg_hh)
#define pylith_meshio_testcellfilteravg_hh

// Include directives ---------------------------------------------------
#include <cppunit/extensions/HelperMacros.h>

// Forward declarations -------------------------------------------------
namespace pylith {
  namespace meshio {
    class TestCellFilterAvg;
  } // meshio
} // pylith

// TestCellFilterAvg ----------------------------------------------------
class pylith::meshio::TestCellFilterAvg : public CppUnit::TestFixture
{ // class TestCellFilterAvg

  // CPPUNIT TEST SUITE /////////////////////////////////////////////////
  CPPUNIT_TEST_SUITE( TestCellFilterAvg );

  CPPUNIT_TEST( testConstructor );
  CPPUNIT_TEST( testFilterMesh );
  CPPUNIT_TEST( testFilterSubMesh );

  CPPUNIT_TEST_SUITE_END();

  // PUBLIC METHODS /////////////////////////////////////////////////////
public :

  /// Test constructor
  void testConstructor(void);

  /// Test filter() w/mesh.
  void testFilterMesh(void);

  /// Test filter() w/submesh.
  void testFilterSubMesh(void);

}; // class TestCellFilterAvg

#endif // pylith_meshio_testcellfilteravg_hh


// End of file 
