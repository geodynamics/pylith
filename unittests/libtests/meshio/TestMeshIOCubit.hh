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
 * @file unittests/libtests/meshio/TestMeshIOCubit.hh
 *
 * @brief C++ TestMeshIOCubit object
 *
 * C++ unit testing for MeshIOCubit.
 */

#if !defined(pylith_meshio_testmeshiocubit_hh)
#define pylith_meshio_testmeshiocubit_hh

#include "TestMeshIO.hh"

/// Namespace for pylith package
namespace pylith {
  namespace meshio {
    class TestMeshIOCubit;
    class MeshData;
  } // meshio
} // pylith

/// C++ unit testing for Quadrature1D
class pylith::meshio::TestMeshIOCubit : public TestMeshIO
{ // class TestMeshIOCubit

  // CPPUNIT TEST SUITE /////////////////////////////////////////////////
  CPPUNIT_TEST_SUITE( TestMeshIOCubit );
  CPPUNIT_TEST( testConstructor );
  CPPUNIT_TEST( testDebug );
  CPPUNIT_TEST( testInterpolate );
  CPPUNIT_TEST( testFilename );
  CPPUNIT_TEST( testReadTet );
  CPPUNIT_TEST( testReadHex );
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

  /// Test read() for mesh with tetrahedral cells.
  void testReadTet(void);

  /// Test read() for mesh with hexahedral cells.
  void testReadHex(void);

  // PRIVATE METHODS ////////////////////////////////////////////////////
private :

  /** Perform read() and then check values.
   *
   * @param data Mesh data
   * @param filename Name of mesh file to read
   */
  void _testRead(const MeshData& data,
		 const char* filename);

}; // class TestMeshIOCubit

#endif // pylith_meshio_testmeshiocubit_hh


// End of file 
