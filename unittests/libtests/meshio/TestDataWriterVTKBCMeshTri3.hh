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
 * @file unittests/libtests/meshio/TestDataWriterVTKBCMeshTri3.hh
 *
 * @brief C++ TestDataWriterVTKBCMeshTri3 object
 *
 * C++ unit testing for DataWriterVTKBCMeshTri3.
 */

#if !defined(pylith_meshio_testdatawritervtkbcmeshtri3_hh)
#define pylith_meshio_testdatawritervtkbcmeshtri3_hh

#include "TestDataWriterVTKBCMesh.hh"

/// Namespace for pylith package
namespace pylith {
  namespace meshio {
    class TestDataWriterVTKBCMeshTri3;
  } // meshio
} // pylith

/// C++ unit testing for DataWriterVTK
class pylith::meshio::TestDataWriterVTKBCMeshTri3 : public TestDataWriterVTKBCMesh
{ // class TestDataWriterVTKBCMeshTri3

  // CPPUNIT TEST SUITE /////////////////////////////////////////////////
  CPPUNIT_TEST_SUITE( TestDataWriterVTKBCMeshTri3 );

  CPPUNIT_TEST( testTimeStep );
  CPPUNIT_TEST( testWriteVertexField );
  CPPUNIT_TEST( testWriteCellField );

  CPPUNIT_TEST_SUITE_END();

  // PUBLIC METHODS /////////////////////////////////////////////////////
public :

  /// Setup testing data.
  void setUp(void);

}; // class TestDataWriterVTKBCMeshTri3

#endif // pylith_meshio_testdatawritervtkbcmeshtri3_hh


// End of file 
