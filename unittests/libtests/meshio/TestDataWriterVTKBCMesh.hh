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
 * @file unittests/libtests/meshio/TestDataWriterVTKBCMesh.hh
 *
 * @brief C++ TestDataWriterVTKBCMesh object
 *
 * C++ unit testing for DataWriterVTKBCMesh.
 */

#if !defined(pylith_meshio_testdatawritervtkbcmesh_hh)
#define pylith_meshio_testdatawritervtkbcmesh_hh

#include "TestDataWriterVTK.hh"

/// Namespace for pylith package
namespace pylith {
  namespace meshio {
    class TestDataWriterVTKBCMesh;

    class DataWriterVTKDataMesh;
  } // meshio
} // pylith

/// C++ unit testing for DataWriterVTK
class pylith::meshio::TestDataWriterVTKBCMesh : public TestDataWriterVTK
{ // class TestDataWriterVTKBCMesh

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

  DataWriterVTKDataMesh* _dataMesh; ///< Data for testing.
  ALE::Obj<Mesh> _meshDomain; ///< Mesh for domain.

}; // class TestDataWriterVTKBCMesh

#endif // pylith_meshio_testdatawritervtkbcmesh_hh


// End of file 
