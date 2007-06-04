// -*- C++ -*-
//
// ----------------------------------------------------------------------
//
//                           Brad T. Aagaard
//                        U.S. Geological Survey
//
// {LicenseText}
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

#include "TestMeshIO.hh"

/// Namespace for pylith package
namespace pylith {
  namespace meshio {
    class TestMeshIOLagrit;
    class MeshData;
  } // meshio
} // pylith

/// C++ unit testing for Quadrature1D
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
   */
  void _testRead(const MeshData& data,
		 const char* filenameGmv,
		 const char* filenamePset);

}; // class TestMeshIOLagrit

#endif // pylith_meshio_testmeshiolagrit_hh


// End of file 
