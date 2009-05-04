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

#include "pylith/topology/topologyfwd.hh" // USES Mesh, SubMesh, Field

/// Namespace for pylith package
namespace pylith {
  namespace meshio {
    class TestDataWriterVTKBCMesh;
  } // meshio
} // pylith

/// C++ unit testing for DataWriterVTK
class pylith::meshio::TestDataWriterVTKBCMesh : public TestDataWriterVTK
{ // class TestDataWriterVTKBCMesh

  // CPPUNIT TEST SUITE /////////////////////////////////////////////////
  CPPUNIT_TEST_SUITE( TestDataWriterVTKBCMesh );

  CPPUNIT_TEST( testConstructor );

  CPPUNIT_TEST_SUITE_END();

  // PUBLIC METHODS /////////////////////////////////////////////////////
public :

  /// Setup testing data.
  void setUp(void);

  /// Tear down testing data.
  void tearDown(void);

  /// Test constructor
  void testConstructor(void);

  /// Test openTimeStep() and closeTimeStep()
  void testTimeStep(void);

  /// Test writeVertexField.
  void testWriteVertexField(void);

  /// Test writeCellField.
  void testWriteCellField(void);

  // PROTECTED MEMBERS //////////////////////////////////////////////////
protected :

  /// Initialize mesh.
  void _initialize(void);

  // PROTECTED MEMBERS //////////////////////////////////////////////////
protected :

  topology::Mesh* _mesh; ///< Mesh for domain
  topology::SubMesh* _submesh; ///< Mesh for subdomain.
  bool _flipFault; ///< If true, flip fault orientation.

  // PRIVATE MEMBERS ////////////////////////////////////////////////////
private :

  /** Create vertex fields.
   *
   * @param fields Vertex fields.
   */
  void
  _createVertexFields(topology::Fields<topology::Field<topology::SubMesh> >* fields) const;

  /** Create cell fields.
   *
   * @param fields Cell fields.
   */
  void
  _createCellFields(topology::Fields<topology::Field<topology::SubMesh> >* fields) const;

}; // class TestDataWriterVTKBCMesh

#endif // pylith_meshio_testdatawritervtkbcmesh_hh


// End of file 
