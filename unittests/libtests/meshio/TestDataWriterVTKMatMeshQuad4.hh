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
 * @file unittests/libtests/meshio/TestDataWriterVTKMatMeshQuad4.hh
 *
 * @brief C++ TestDataWriterVTKMatMeshQuad4 object
 *
 * C++ unit testing for DataWriterVTKMatMeshQuad4.
 */

#if !defined(pylith_meshio_testdatawritervtksubmeshquad4_hh)
#define pylith_meshio_testdatawritervtksubmeshquad4_hh

#include "TestDataWriterVTKMesh.hh"

/// Namespace for pylith package
namespace pylith {
  namespace meshio {
    class TestDataWriterVTKMatMeshQuad4;
  } // meshio
} // pylith

/// C++ unit testing for DataWriterVTK
class pylith::meshio::TestDataWriterVTKMatMeshQuad4 : public TestDataWriterVTKMesh
{ // class TestDataWriterVTKMatMeshQuad4

  // CPPUNIT TEST SUITE /////////////////////////////////////////////////
  CPPUNIT_TEST_SUITE( TestDataWriterVTKMatMeshQuad4 );

  CPPUNIT_TEST( testTimeStep );
  CPPUNIT_TEST( testWriteVertexField );
  CPPUNIT_TEST( testWriteCellField );

  CPPUNIT_TEST_SUITE_END();

  // PUBLIC METHODS /////////////////////////////////////////////////////
public :

  /// Setup testing data.
  void setUp(void);

}; // class TestDataWriterVTKMatMeshQuad4

#endif // pylith_meshio_testdatawritervtksubmeshquad4_hh


// End of file 
