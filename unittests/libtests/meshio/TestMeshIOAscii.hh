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
 * @file unittests/libtests/meshio/TestMeshIOAscii.hh
 *
 * @brief C++ TestMeshIOAscii object
 *
 * C++ unit testing for MeshIOAscii.
 */

#if !defined(pylith_meshio_testmeshioascii_hh)
#define pylith_meshio_testmeshioascii_hh

// Include directives ---------------------------------------------------
#include "TestMeshIO.hh"

// Forward declarations -------------------------------------------------
namespace pylith {
  namespace meshio {
    class TestMeshIOAscii;
    class MeshData;
  } // meshio
} // pylith

// TestMeshIOAscii ------------------------------------------------------
class pylith::meshio::TestMeshIOAscii : public TestMeshIO
{ // class TestMeshIOAscii

  // CPPUNIT TEST SUITE /////////////////////////////////////////////////
  CPPUNIT_TEST_SUITE( TestMeshIOAscii );

  CPPUNIT_TEST( testConstructor );
  CPPUNIT_TEST( testDebug );
  CPPUNIT_TEST( testInterpolate );
  CPPUNIT_TEST( testFilename );
  CPPUNIT_TEST( testWriteRead1D );
  CPPUNIT_TEST( testWriteRead1Din2D );
  CPPUNIT_TEST( testWriteRead1Din3D );
  CPPUNIT_TEST( testWriteRead2D );
  CPPUNIT_TEST( testWriteRead2Din3D );
  CPPUNIT_TEST( testWriteRead3D );
  CPPUNIT_TEST( testRead3DIndexOne );
  CPPUNIT_TEST( testReadComments );

  CPPUNIT_TEST_SUITE_END();

  // PUBLIC METHODS /////////////////////////////////////////////////////
public :

  /// Test constructor
  void testConstructor(void);

  /// Test debug()
  void testDebug(void);

  /// Test interpolate()
  void testInterpolate(void);

  /// Test filename()
  void testFilename(void);

  /// Test write() and read() for 1D mesh in 1D space.
  void testWriteRead1D(void);

  /// Test write() and read() for 1D mesh in 2D space.
  void testWriteRead1Din2D(void);

  /// Test write() and read() for 1D mesh in 3D space.
  void testWriteRead1Din3D(void);

  /// Test write() and read() for 2D mesh in 2D space.
  void testWriteRead2D(void);

  /// Test write() and read() for 2D mesh in 3D space.
  void testWriteRead2Din3D(void);

  /// Test write() and read() for 3D mesh in 3D space.
  void testWriteRead3D(void);

  /// Test read() for 3D mesh with one based indexing.
  void testRead3DIndexOne(void);

  /// Test and read() for 2D mesh in 2D space with comments.
  void testReadComments(void);

  // PRIVATE METHODS ////////////////////////////////////////////////////
private :

  /** Build mesh, perform write() and read(), and then check values.
   *
   * @param data Mesh data
   * @param filename Name of mesh file to write/read
   */
  void _testWriteRead(const MeshData& data,
		      const char* filename);

  /** Read mesh and then check values.
   *
   * @param data Mesh data
   * @param filename Name of mesh file to read
   */
  void _testRead(const MeshData& data,
		 const char* filename);

}; // class TestMeshIOAscii

#endif // pylith_meshio_testmeshioascii_hh

// End of file 
