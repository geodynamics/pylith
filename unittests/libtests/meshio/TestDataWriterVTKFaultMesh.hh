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
 * @file unittests/libtests/meshio/TestDataWriterVTKFaultMesh.hh
 *
 * @brief C++ TestDataWriterVTKFaultMesh object
 *
 * C++ unit testing for DataWriterVTKFaultMesh.
 */

#if !defined(pylith_meshio_testdatawritervtkfaultmesh_hh)
#define pylith_meshio_testdatawritervtkfaultmesh_hh

#include "TestDataWriterVTK.hh"

#include "pylith/topology/topologyfwd.hh" // USES Mesh, SubMesh, Field

/// Namespace for pylith package
namespace pylith {
  namespace meshio {
    class TestDataWriterVTKFaultMesh;
  } // meshio
} // pylith

/// C++ unit testing for DataWriterVTK
class pylith::meshio::TestDataWriterVTKFaultMesh : public TestDataWriterVTK
{ // class TestDataWriterVTKFaultMesh

  // CPPUNIT TEST SUITE /////////////////////////////////////////////////
  CPPUNIT_TEST_SUITE( TestDataWriterVTKFaultMesh );

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
  topology::SubMesh* _faultMesh; ///< Fault mesh.
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

}; // class TestDataWriterVTKFaultMesh

#endif // pylith_meshio_testdatawritervtkfaultmesh_hh


// End of file 
