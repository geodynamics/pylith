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
 * @file unittests/libtests/meshio/TestDataWriterVTKSubMeshHex8.hh
 *
 * @brief C++ TestDataWriterVTKSubMeshHex8 object
 *
 * C++ unit testing for DataWriterVTKSubMeshHex8.
 */

#if !defined(pylith_meshio_testdatawritervtksubmeshhex8_hh)
#define pylith_meshio_testdatawritervtksubmeshhex8_hh

#include "TestDataWriterVTKMesh.hh"

/// Namespace for pylith package
namespace pylith {
  namespace meshio {
    class TestDataWriterVTKSubMeshHex8;
  } // meshio
} // pylith

/// C++ unit testing for DataWriterVTK
class pylith::meshio::TestDataWriterVTKSubMeshHex8 : public TestDataWriterVTKMesh
{ // class TestDataWriterVTKSubMeshHex8

  // CPPUNIT TEST SUITE /////////////////////////////////////////////////
  CPPUNIT_TEST_SUITE( TestDataWriterVTKSubMeshHex8 );

  CPPUNIT_TEST( testTimeStep );
  CPPUNIT_TEST( testWriteVertexField );
  CPPUNIT_TEST( testWriteCellField );

  CPPUNIT_TEST_SUITE_END();

  // PUBLIC METHODS /////////////////////////////////////////////////////
public :

  /// Setup testing data.
  void setUp(void);

}; // class TestDataWriterVTKSubMeshHex8

#endif // pylith_meshio_testdatawritervtksubmeshhex8_hh


// End of file 
