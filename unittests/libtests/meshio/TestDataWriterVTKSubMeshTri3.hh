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
 * @file unittests/libtests/meshio/TestDataWriterVTKSubMeshTri3.hh
 *
 * @brief C++ TestDataWriterVTKSubMeshTri3 object
 *
 * C++ unit testing for DataWriterVTKSubMeshTri3.
 */

#if !defined(pylith_meshio_testdatawritervtksubmeshtri3_hh)
#define pylith_meshio_testdatawritervtksubmeshtri3_hh

#include "TestDataWriterVTKSubMesh.hh"

/// Namespace for pylith package
namespace pylith {
  namespace meshio {
    class TestDataWriterVTKSubMeshTri3;
  } // meshio
} // pylith

/// C++ unit testing for DataWriterVTK
class pylith::meshio::TestDataWriterVTKSubMeshTri3 : public TestDataWriterVTKSubMesh
{ // class TestDataWriterVTKSubMeshTri3

  // CPPUNIT TEST SUITE /////////////////////////////////////////////////
  CPPUNIT_TEST_SUITE( TestDataWriterVTKSubMeshTri3 );

  CPPUNIT_TEST( testTimeStep );
  CPPUNIT_TEST( testWriteVertexField );
  CPPUNIT_TEST( testWriteCellField );

  CPPUNIT_TEST_SUITE_END();

  // PUBLIC METHODS /////////////////////////////////////////////////////
public :

  /// Setup testing data.
  void setUp(void);

}; // class TestDataWriterVTKSubMeshTri3

#endif // pylith_meshio_testdatawritervtksubmeshtri3_hh


// End of file 
