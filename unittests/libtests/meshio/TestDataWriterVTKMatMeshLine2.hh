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
 * @file unittests/libtests/meshio/TestDataWriterVTKMatMeshLine2.hh
 *
 * @brief C++ TestDataWriterVTKMatMeshLine2 object
 *
 * C++ unit testing for DataWriterVTKMatMeshLine2.
 */

#if !defined(pylith_meshio_testdatawritervtksubmeshline2_hh)
#define pylith_meshio_testdatawritervtksubmeshline2_hh

#include "TestDataWriterVTKMesh.hh"

/// Namespace for pylith package
namespace pylith {
  namespace meshio {
    class TestDataWriterVTKMatMeshLine2;
  } // meshio
} // pylith

/// C++ unit testing for DataWriterVTK
class pylith::meshio::TestDataWriterVTKMatMeshLine2 : public TestDataWriterVTKMesh
{ // class TestDataWriterVTKMatMeshLine2

  // CPPUNIT TEST SUITE /////////////////////////////////////////////////
  CPPUNIT_TEST_SUITE( TestDataWriterVTKMatMeshLine2 );

  CPPUNIT_TEST( testTimeStep );
  CPPUNIT_TEST( testWriteVertexField );
  CPPUNIT_TEST( testWriteCellField );

  CPPUNIT_TEST_SUITE_END();

  // PUBLIC METHODS /////////////////////////////////////////////////////
public :

  /// Setup testing data.
  void setUp(void);

}; // class TestDataWriterVTKMatMeshLine2

#endif // pylith_meshio_testdatawritervtksubmeshline2_hh


// End of file 
