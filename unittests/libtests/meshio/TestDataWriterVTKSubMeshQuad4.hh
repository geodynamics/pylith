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
 * @file unittests/libtests/meshio/TestDataWriterVTKSubMeshQuad4.hh
 *
 * @brief C++ TestDataWriterVTKSubMeshQuad4 object
 *
 * C++ unit testing for DataWriterVTKSubMeshQuad4.
 */

#if !defined(pylith_meshio_testdatawritervtksubmeshquad4_hh)
#define pylith_meshio_testdatawritervtksubmeshquad4_hh

#include "TestDataWriterVTKSubMesh.hh"

/// Namespace for pylith package
namespace pylith {
  namespace meshio {
    class TestDataWriterVTKSubMeshQuad4;
  } // meshio
} // pylith

/// C++ unit testing for DataWriterVTK
class pylith::meshio::TestDataWriterVTKSubMeshQuad4 : public TestDataWriterVTKSubMesh
{ // class TestDataWriterVTKSubMeshQuad4

  // CPPUNIT TEST SUITE /////////////////////////////////////////////////
  CPPUNIT_TEST_SUITE( TestDataWriterVTKSubMeshQuad4 );

  CPPUNIT_TEST( testTimeStep );
  CPPUNIT_TEST( testWriteVertexField );
  CPPUNIT_TEST( testWriteCellField );

  CPPUNIT_TEST_SUITE_END();

  // PUBLIC METHODS /////////////////////////////////////////////////////
public :

  /// Setup testing data.
  void setUp(void);

}; // class TestDataWriterVTKSubMeshQuad4

#endif // pylith_meshio_testdatawritervtksubmeshquad4_hh


// End of file 
