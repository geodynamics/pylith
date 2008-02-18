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
 * @file unittests/libtests/meshio/TestDataWriterVTKMeshQuad4.hh
 *
 * @brief C++ TestDataWriterVTKMeshQuad4 object
 *
 * C++ unit testing for DataWriterVTKMeshQuad4.
 */

#if !defined(pylith_meshio_testdatawritervtkmeshquad4_hh)
#define pylith_meshio_testdatawritervtkmeshquad4_hh

#include "TestDataWriterVTKMesh.hh"

/// Namespace for pylith package
namespace pylith {
  namespace meshio {
    class TestDataWriterVTKMeshQuad4;
  } // meshio
} // pylith

/// C++ unit testing for DataWriterVTK
class pylith::meshio::TestDataWriterVTKMeshQuad4 : public TestDataWriterVTKMesh
{ // class TestDataWriterVTKMeshQuad4

  // CPPUNIT TEST SUITE /////////////////////////////////////////////////
  CPPUNIT_TEST_SUITE( TestDataWriterVTKMeshQuad4 );

  CPPUNIT_TEST( testTimeStep );
  CPPUNIT_TEST( testWriteVertexField );
  CPPUNIT_TEST( testWriteCellField );

  CPPUNIT_TEST_SUITE_END();

  // PUBLIC METHODS /////////////////////////////////////////////////////
public :

  /// Setup testing data.
  void setUp(void);

}; // class TestDataWriterVTKMeshQuad4

#endif // pylith_meshio_testdatawritervtkmeshquad4_hh


// End of file 
