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
 * @file unittests/libtests/meshio/TestDataWriterVTKFaultMeshTri3.hh
 *
 * @brief C++ TestDataWriterVTKFaultMeshTri3 object
 *
 * C++ unit testing for DataWriterVTKFaultMeshTri3.
 */

#if !defined(pylith_meshio_testdatawritervtkfaultmeshtri3_hh)
#define pylith_meshio_testdatawritervtkfaultmeshtri3_hh

#include "TestDataWriterVTKFaultMesh.hh"

/// Namespace for pylith package
namespace pylith {
  namespace meshio {
    class TestDataWriterVTKFaultMeshTri3;
  } // meshio
} // pylith

/// C++ unit testing for DataWriterVTK
class pylith::meshio::TestDataWriterVTKFaultMeshTri3 : public TestDataWriterVTKFaultMesh
{ // class TestDataWriterVTKFaultMeshTri3

  // CPPUNIT TEST SUITE /////////////////////////////////////////////////
  CPPUNIT_TEST_SUITE( TestDataWriterVTKFaultMeshTri3 );

  CPPUNIT_TEST( testTimeStep );
  CPPUNIT_TEST( testWriteVertexField );
  CPPUNIT_TEST( testWriteCellField );

  CPPUNIT_TEST_SUITE_END();

  // PUBLIC METHODS /////////////////////////////////////////////////////
public :

  /// Setup testing data.
  void setUp(void);

}; // class TestDataWriterVTKFaultMeshTri3

#endif // pylith_meshio_testdatawritervtkfaultmeshtri3_hh


// End of file 
