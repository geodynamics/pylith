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
 * @file unittests/libtests/meshio/TestMeshIOLagrit.hh
 *
 * @brief C++ TestMeshIOLagrit object
 *
 * C++ unit testing for MeshIOLagrit.
 */

#if !defined(pylith_meshio_testmeshiolagrit_hh)
#define pylith_meshio_testmeshiolagrit_hh

// Include directives ---------------------------------------------------
#include "TestMeshIO.hh"

// Forward declarations -------------------------------------------------
namespace pylith {
  namespace meshio {
    class TestMeshIOLagrit;
    class MeshData;
  } // meshio
} // pylith

// TestMeshIOLagrit -----------------------------------------------------
class pylith::meshio::TestMeshIOLagrit : public TestMeshIO
{ // class TestMeshIOLagrit

  // CPPUNIT TEST SUITE /////////////////////////////////////////////////
  CPPUNIT_TEST_SUITE( TestMeshIOLagrit );
  CPPUNIT_TEST( testConstructor );
  CPPUNIT_TEST( testDebug );
  CPPUNIT_TEST( testInterpolate );
  CPPUNIT_TEST( testFilename );
  CPPUNIT_TEST( testReadTetAscii );
  CPPUNIT_TEST( testReadTetBinary );
  CPPUNIT_TEST( testReadTetBinary32on64 );
  CPPUNIT_TEST( testReadTetBinary64 );
  CPPUNIT_TEST( testOrientAsciiTet );
  CPPUNIT_TEST( testOrientBinaryTet );
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

  /// Test read() for mesh with ASCII files.
  void testReadTetAscii(void);

  /// Test read() for mesh with binary files.
  void testReadTetBinary(void);

  /// Test read() for mesh with binary files for 32-bit LaGrit built 
  /// on 64-bit platform.
  void testReadTetBinary32on64(void);

  /// Test read() for mesh with binary files for 64-bit LaGriT.
  void testReadTetBinary64(void);

  /// Test _orientCellsAscii with tet cells.
  void testOrientAsciiTet(void);

  /// Test _orientCellsBinary with tet cells.
  void testOrientBinaryTet(void);

  // PRIVATE METHODS ////////////////////////////////////////////////////
private :

  /** Perform read() and then check values.
   *
   * @param data Mesh data
   * @param filenameGmv Name of mesh GMV file to read
   * @param filenamePset Name of mesh Pset file to read
   * @param ioInt32 True if Pset uses 32-bit integers.
   * @param isRecordHeader32Bit True if Fortran record headers are 32-bit.
   */
  void _testRead(const MeshData& data,
		 const char* filenameGmv,
		 const char* filenamePset,
		 const bool ioInt32 =true,
		 const bool isRecordHeader32Bit =true);

}; // class TestMeshIOLagrit

#endif // pylith_meshio_testmeshiolagrit_hh


// End of file 
