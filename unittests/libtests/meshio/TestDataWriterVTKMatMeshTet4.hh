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
 * @file unittests/libtests/meshio/TestDataWriterVTKMatMeshTet4.hh
 *
 * @brief C++ TestDataWriterVTKMatMeshTet4 object
 *
 * C++ unit testing for DataWriterVTKMatMeshTet4.
 */

#if !defined(pylith_meshio_testdatawritervtksubmeshtet4_hh)
#define pylith_meshio_testdatawritervtksubmeshtet4_hh

#include "TestDataWriterVTKMesh.hh"

/// Namespace for pylith package
namespace pylith {
  namespace meshio {
    class TestDataWriterVTKMatMeshTet4;
  } // meshio
} // pylith

/// C++ unit testing for DataWriterVTK
class pylith::meshio::TestDataWriterVTKMatMeshTet4 : public TestDataWriterVTKMesh
{ // class TestDataWriterVTKMatMeshTet4

  // CPPUNIT TEST SUITE /////////////////////////////////////////////////
  CPPUNIT_TEST_SUITE( TestDataWriterVTKMatMeshTet4 );

  CPPUNIT_TEST( testTimeStep );
  CPPUNIT_TEST( testWriteVertexField );
  CPPUNIT_TEST( testWriteCellField );

  CPPUNIT_TEST_SUITE_END();

  // PUBLIC METHODS /////////////////////////////////////////////////////
public :

  /// Setup testing data.
  void setUp(void);

}; // class TestDataWriterVTKMatMeshTet4

#endif // pylith_meshio_testdatawritervtksubmeshtet4_hh


// End of file 
