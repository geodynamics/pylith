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
 * @file unittests/libtests/meshio/TestDataWriterVTKMesh.hh
 *
 * @brief C++ TestDataWriterVTKMesh object
 *
 * C++ unit testing for DataWriterVTKMesh.
 */

#if !defined(pylith_meshio_testdatawritervtkmesh_hh)
#define pylith_meshio_testdatawritervtkmesh_hh

#include "TestDataWriterVTK.hh"

/// Namespace for pylith package
namespace pylith {
  namespace meshio {
    class TestDataWriterVTKMesh;

    class DataWriterVTKDataMesh;
  } // meshio
} // pylith

/// C++ unit testing for DataWriterVTK
class pylith::meshio::TestDataWriterVTKMesh : public TestDataWriterVTK
{ // class TestDataWriterVTKMesh

  // PUBLIC METHODS /////////////////////////////////////////////////////
public :

  /// Setup testing data.
  void setUp(void);

  /// Tear down testing data.
  void tearDown(void);

  // PROTECTED MEMBERS //////////////////////////////////////////////////
protected :

  /// Initialize mesh.
  void _initialize(void);

  // PROTECTED MEMBERS //////////////////////////////////////////////////
protected :

  DataWriterVTKDataMesh* _dataMesh; ///< Data for testing

}; // class TestDataWriterVTKMesh

#endif // pylith_meshio_testdatawritervtkmesh_hh


// End of file 
