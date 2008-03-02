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
 * @file unittests/libtests/meshio/TestDataWriterVTKBCMeshQuad4.hh
 *
 * @brief C++ TestDataWriterVTKBCMeshQuad4 object
 *
 * C++ unit testing for DataWriterVTKBCMeshQuad4.
 */

#if !defined(pylith_meshio_testdatawritervtkbcmeshquad4_hh)
#define pylith_meshio_testdatawritervtkbcmeshquad4_hh

#include "TestDataWriterVTKBCMesh.hh"

/// Namespace for pylith package
namespace pylith {
  namespace meshio {
    class TestDataWriterVTKBCMeshQuad4;
  } // meshio
} // pylith

/// C++ unit testing for DataWriterVTK
class pylith::meshio::TestDataWriterVTKBCMeshQuad4 : public TestDataWriterVTKBCMesh
{ // class TestDataWriterVTKBCMeshQuad4

  // CPPUNIT TEST SUITE /////////////////////////////////////////////////
  CPPUNIT_TEST_SUITE( TestDataWriterVTKBCMeshQuad4 );

  CPPUNIT_TEST( testTimeStep );
  CPPUNIT_TEST( testWriteVertexField );
  CPPUNIT_TEST( testWriteCellField );

  CPPUNIT_TEST_SUITE_END();

  // PUBLIC METHODS /////////////////////////////////////////////////////
public :

  /// Setup testing data.
  void setUp(void);

}; // class TestDataWriterVTKBCMeshQuad4

#endif // pylith_meshio_testdatawritervtkbcmeshquad4_hh


// End of file 
