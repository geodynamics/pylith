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
 * @file unittests/libtests/meshio/TestMeshIOAscii.hh
 *
 * @brief C++ TestMeshIOAscii object
 *
 * C++ unit testing for MeshIOAscii.
 */

#if !defined(pylith_meshio_testmeshioascii_hh)
#define pylith_meshio_testmeshioascii_hh

#include "TestMeshIO.hh"

/// Namespace for pylith package
namespace pylith {
  namespace meshio {
    class TestMeshIOAscii;
    class MeshData;
  } // meshio
} // pylith

/// C++ unit testing for Quadrature1D
class pylith::meshio::TestMeshIOAscii : public TestMeshIO
{ // class TestMeshIOAscii

  // CPPUNIT TEST SUITE /////////////////////////////////////////////////
  CPPUNIT_TEST_SUITE( TestMeshIOAscii );
  CPPUNIT_TEST( testConstructor );
  CPPUNIT_TEST( testFilename );
  CPPUNIT_TEST( testWriteRead1D );
  CPPUNIT_TEST( testWriteRead1Din2D );
  CPPUNIT_TEST( testWriteRead1Din3D );
  CPPUNIT_TEST( testWriteRead2D );
  CPPUNIT_TEST( testWriteRead2Din3D );
  CPPUNIT_TEST( testWriteRead3D );
  CPPUNIT_TEST_SUITE_END();

  // PUBLIC METHODS /////////////////////////////////////////////////////
public :

  /// Test constructor
  void testConstructor(void);

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

  // PRIVATE METHODS ////////////////////////////////////////////////////
private :

  /** Build mesh, perform write() and read(), and then check values.
   *
   * @param data Mesh data
   * @param filename Name of mesh file to write/read
   */
  void _testWriteRead(const MeshData& data,
		      const char* filename);

}; // class TestMeshIOAscii

#endif // pylith_meshio_testmeshioascii_hh

// End of file 
