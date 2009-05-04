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
 * @file unittests/libtests/meshio/TestDataWriterVTKSubMeshLine2.hh
 *
 * @brief C++ TestDataWriterVTKSubMeshLine2 object
 *
 * C++ unit testing for DataWriterVTKSubMeshLine2.
 */

#if !defined(pylith_meshio_testdatawritervtksubmeshline2_hh)
#define pylith_meshio_testdatawritervtksubmeshline2_hh

#include "TestDataWriterVTKSubMesh.hh"

/// Namespace for pylith package
namespace pylith {
  namespace meshio {
    class TestDataWriterVTKSubMeshLine2;
  } // meshio
} // pylith

/// C++ unit testing for DataWriterVTK
class pylith::meshio::TestDataWriterVTKSubMeshLine2 : public TestDataWriterVTKSubMesh
{ // class TestDataWriterVTKSubMeshLine2

  // CPPUNIT TEST SUITE /////////////////////////////////////////////////
  CPPUNIT_TEST_SUITE( TestDataWriterVTKSubMeshLine2 );

  CPPUNIT_TEST( testTimeStep );
  CPPUNIT_TEST( testWriteVertexField );
  CPPUNIT_TEST( testWriteCellField );

  CPPUNIT_TEST_SUITE_END();

  // PUBLIC METHODS /////////////////////////////////////////////////////
public :

  /// Setup testing data.
  void setUp(void);

}; // class TestDataWriterVTKSubMeshLine2

#endif // pylith_meshio_testdatawritervtksubmeshline2_hh


// End of file 
