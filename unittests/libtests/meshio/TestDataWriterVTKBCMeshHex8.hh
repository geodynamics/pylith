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
 * @file unittests/libtests/meshio/TestDataWriterVTKBCMeshHex8.hh
 *
 * @brief C++ TestDataWriterVTKBCMeshHex8 object
 *
 * C++ unit testing for DataWriterVTKBCMeshHex8.
 */

#if !defined(pylith_meshio_testdatawritervtkbcmeshhex8_hh)
#define pylith_meshio_testdatawritervtkbcmeshhex8_hh

#include "TestDataWriterVTKBCMesh.hh"

/// Namespace for pylith package
namespace pylith {
  namespace meshio {
    class TestDataWriterVTKBCMeshHex8;
  } // meshio
} // pylith

/// C++ unit testing for DataWriterVTK
class pylith::meshio::TestDataWriterVTKBCMeshHex8 : public TestDataWriterVTKBCMesh
{ // class TestDataWriterVTKBCMeshHex8

  // CPPUNIT TEST SUITE /////////////////////////////////////////////////
  CPPUNIT_TEST_SUITE( TestDataWriterVTKBCMeshHex8 );

  CPPUNIT_TEST( testTimeStep );
  CPPUNIT_TEST( testWriteVertexField );
  CPPUNIT_TEST( testWriteCellField );

  CPPUNIT_TEST_SUITE_END();

  // PUBLIC METHODS /////////////////////////////////////////////////////
public :

  /// Setup testing data.
  void setUp(void);

}; // class TestDataWriterVTKBCMeshHex8

#endif // pylith_meshio_testdatawritervtkbcmeshhex8_hh


// End of file 
