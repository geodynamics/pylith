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
 * @file unittests/libtests/meshio/TestDataWriterVTKFaultMeshQuad4.hh
 *
 * @brief C++ TestDataWriterVTKFaultMeshQuad4 object
 *
 * C++ unit testing for DataWriterVTKFaultMeshQuad4.
 */

#if !defined(pylith_meshio_testdatawritervtkfaultmeshquad4_hh)
#define pylith_meshio_testdatawritervtkfaultmeshquad4_hh

#include "TestDataWriterVTKFaultMesh.hh"

/// Namespace for pylith package
namespace pylith {
  namespace meshio {
    class TestDataWriterVTKFaultMeshQuad4;
  } // meshio
} // pylith

/// C++ unit testing for DataWriterVTK
class pylith::meshio::TestDataWriterVTKFaultMeshQuad4 : public TestDataWriterVTKFaultMesh
{ // class TestDataWriterVTKFaultMeshQuad4

  // CPPUNIT TEST SUITE /////////////////////////////////////////////////
  CPPUNIT_TEST_SUITE( TestDataWriterVTKFaultMeshQuad4 );

  CPPUNIT_TEST( testTimeStep );
  CPPUNIT_TEST( testWriteVertexField );
  CPPUNIT_TEST( testWriteCellField );

  CPPUNIT_TEST_SUITE_END();

  // PUBLIC METHODS /////////////////////////////////////////////////////
public :

  /// Setup testing data.
  void setUp(void);

}; // class TestDataWriterVTKFaultMeshQuad4

#endif // pylith_meshio_testdatawritervtkfaultmeshquad4_hh


// End of file 
