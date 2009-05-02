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
 * @file unittests/libtests/meshio/TestDataWriterVTKMatMeshTri3.hh
 *
 * @brief C++ TestDataWriterVTKMatMeshTri3 object
 *
 * C++ unit testing for DataWriterVTKMatMeshTri3.
 */

#if !defined(pylith_meshio_testdatawritervtksubmeshtri3_hh)
#define pylith_meshio_testdatawritervtksubmeshtri3_hh

#include "TestDataWriterVTKMesh.hh"

/// Namespace for pylith package
namespace pylith {
  namespace meshio {
    class TestDataWriterVTKMatMeshTri3;
  } // meshio
} // pylith

/// C++ unit testing for DataWriterVTK
class pylith::meshio::TestDataWriterVTKMatMeshTri3 : public TestDataWriterVTKMesh
{ // class TestDataWriterVTKMatMeshTri3

  // CPPUNIT TEST SUITE /////////////////////////////////////////////////
  CPPUNIT_TEST_SUITE( TestDataWriterVTKMatMeshTri3 );

  CPPUNIT_TEST( testTimeStep );
  CPPUNIT_TEST( testWriteVertexField );
  CPPUNIT_TEST( testWriteCellField );

  CPPUNIT_TEST_SUITE_END();

  // PUBLIC METHODS /////////////////////////////////////////////////////
public :

  /// Setup testing data.
  void setUp(void);

}; // class TestDataWriterVTKMatMeshTri3

#endif // pylith_meshio_testdatawritervtksubmeshtri3_hh


// End of file 
