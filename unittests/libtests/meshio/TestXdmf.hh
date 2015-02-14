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
  CPPUNIT_TEST( testWriteTri3Vertex );
  CPPUNIT_TEST( testWriteTri3Cell );
  CPPUNIT_TEST( testWriteQuad4Vertex );
  CPPUNIT_TEST( testWriteQuad4Cell );
  CPPUNIT_TEST( testWriteTet4Vertex );
  CPPUNIT_TEST( testWriteTet4Cell );
  CPPUNIT_TEST( testWriteHex8Vertex );
  CPPUNIT_TEST( testWriteHex8Cell );

  CPPUNIT_TEST_SUITE_END();

  // PUBLIC METHODS /////////////////////////////////////////////////////
public :

  /// Test constructor.
  void testConstructor(void);

  /// Test write() with tri3 mesh and vertex data.
  void testWriteTri3Vertex(void);

  /// Test write() with tri3 mesh and cell data.
  void testWriteTri3Cell(void);

  /// Test write() with quad4 mesh and vertex data.
  void testWriteQuad4Vertex(void);

  /// Test write() with quad4 mesh and cell data.
  void testWriteQuad4Cell(void);

  /// Test write() with tet4 mesh and vertex data.
  void testWriteTet4Vertex(void);

  /// Test write() with tet4 mesh and cell data.
  void testWriteTet4Cell(void);

  /// Test write() with hex8 mesh and vertex data.
  void testWriteHex8Vertex(void);

  /// Test write() with hex8 mesh and cell data.
  void testWriteHex8Cell(void);

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
