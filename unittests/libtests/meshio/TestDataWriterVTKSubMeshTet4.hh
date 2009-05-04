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
 * @file unittests/libtests/meshio/TestDataWriterVTKSubMeshTet4.hh
 *
 * @brief C++ TestDataWriterVTKSubMeshTet4 object
 *
 * C++ unit testing for DataWriterVTKSubMeshTet4.
 */

#if !defined(pylith_meshio_testdatawritervtksubmeshtet4_hh)
#define pylith_meshio_testdatawritervtksubmeshtet4_hh

#include "TestDataWriterVTKSubMesh.hh"

/// Namespace for pylith package
namespace pylith {
  namespace meshio {
    class TestDataWriterVTKSubMeshTet4;
  } // meshio
} // pylith

/// C++ unit testing for DataWriterVTK
class pylith::meshio::TestDataWriterVTKSubMeshTet4 : public TestDataWriterVTKSubMesh
{ // class TestDataWriterVTKSubMeshTet4

  // CPPUNIT TEST SUITE /////////////////////////////////////////////////
  CPPUNIT_TEST_SUITE( TestDataWriterVTKSubMeshTet4 );

  CPPUNIT_TEST( testTimeStep );
  CPPUNIT_TEST( testWriteVertexField );
  CPPUNIT_TEST( testWriteCellField );

  CPPUNIT_TEST_SUITE_END();

  // PUBLIC METHODS /////////////////////////////////////////////////////
public :

  /// Setup testing data.
  void setUp(void);

}; // class TestDataWriterVTKSubMeshTet4

#endif // pylith_meshio_testdatawritervtksubmeshtet4_hh


// End of file 
