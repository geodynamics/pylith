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
 * @file unittests/libtests/meshio/TestDataWriterVTKMeshTri3.hh
 *
 * @brief C++ TestDataWriterVTKMeshTri3 object
 *
 * C++ unit testing for DataWriterVTKMeshTri3.
 */

#if !defined(pylith_meshio_testdatawritervtkmeshtri3_hh)
#define pylith_meshio_testdatawritervtkmeshtri3_hh

#include "TestDataWriterVTKMesh.hh"

/// Namespace for pylith package
namespace pylith {
  namespace meshio {
    class TestDataWriterVTKMeshTri3;
  } // meshio
} // pylith

/// C++ unit testing for DataWriterVTK
class pylith::meshio::TestDataWriterVTKMeshTri3 : public TestDataWriterVTKMesh
{ // class TestDataWriterVTKMeshTri3

  // CPPUNIT TEST SUITE /////////////////////////////////////////////////
  CPPUNIT_TEST_SUITE( TestDataWriterVTKMeshTri3 );

  CPPUNIT_TEST( testTimeStep );
  CPPUNIT_TEST( testWriteVertexField );
  CPPUNIT_TEST( testWriteCellField );

  CPPUNIT_TEST_SUITE_END();

  // PUBLIC METHODS /////////////////////////////////////////////////////
public :

  /// Setup testing data.
  void setUp(void);

}; // class TestDataWriterVTKMeshTri3

#endif // pylith_meshio_testdatawritervtkmeshtri3_hh


// End of file 
