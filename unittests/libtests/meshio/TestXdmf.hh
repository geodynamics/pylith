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
 * @file unittests/libtests/meshio/TestXdmf.hh
 *
 * @brief C++ TestXdmf object
 *
 * C++ unit testing for Xdmf.
 */

#if !defined(pylith_meshio_testxdmf_hh)
#define pylith_meshio_testxdmf_hh

#include <cppunit/extensions/HelperMacros.h>

/// Namespace for pylith package
namespace pylith {
  namespace meshio {
    class TestXdmf;
  } // meshio
} // pylith

/// C++ unit testing for Xdmf
class pylith::meshio::TestXdmf : public CppUnit::TestFixture
{ // class TestXdmf

  // CPPUNIT TEST SUITE /////////////////////////////////////////////////
  CPPUNIT_TEST_SUITE( TestXdmf );

  CPPUNIT_TEST( testConstructor );
  CPPUNIT_TEST( testWrite2DVertex );
  CPPUNIT_TEST( testWrite2DCell );
  CPPUNIT_TEST( testWrite3DVertex );
  CPPUNIT_TEST( testWrite3DCell );

  CPPUNIT_TEST_SUITE_END();

  // PUBLIC METHODS /////////////////////////////////////////////////////
public :

  /// Test constructor.
  void testConstructor(void);

  /// Test write() with 2D mesh and vertex data.
  void testWrite2DVertex(void);

  /// Test write() with 2D mesh and cell data.
  void testWrite2DCell(void);

  /// Test write() with 3D mesh and vertex data.
  void testWrite3DVertex(void);

  /// Test write() with 3D mesh and cell data.
  void testWrite3DCell(void);

  // PRIVATE METHODS ////////////////////////////////////////////////////
private :

  /** Check Xmdf file against archived file.
   *
   * @param filename Name of file to check.
   */
  static
  void _checkFile(const char* filename);
  
}; // class TestXdmf

#endif // pylith_meshio_testxdmf_hh

// End of file 
