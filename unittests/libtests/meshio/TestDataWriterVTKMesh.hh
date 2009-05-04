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

#include "TestDataWriterVTK.hh" // ISA TestDataWriterVTK

#include "pylith/topology/topologyfwd.hh" // USES Mesh, Field

/// Namespace for pylith package
namespace pylith {
  namespace meshio {
    class TestDataWriterVTKMesh;
  } // meshio
} // pylith

/// C++ unit testing for DataWriterVTK
class pylith::meshio::TestDataWriterVTKMesh : public TestDataWriterVTK
{ // class TestDataWriterVTKMesh

  // CPPUNIT TEST SUITE /////////////////////////////////////////////////
  CPPUNIT_TEST_SUITE( TestDataWriterVTKMesh );

  CPPUNIT_TEST( testConstructor );
  CPPUNIT_TEST( testFilename );
  CPPUNIT_TEST( testTimeFormat );
  CPPUNIT_TEST( testTimeConstant );
  CPPUNIT_TEST( testVtkFilename );

  CPPUNIT_TEST_SUITE_END();

  // PUBLIC METHODS /////////////////////////////////////////////////////
public :

  /// Setup testing data.
  void setUp(void);

  /// Tear down testing data.
  void tearDown(void);

  /// Test constructor
  void testConstructor(void);

  /// Test filename()
  void testFilename(void);

  /// Test timeFormat()
  void testTimeFormat(void);

  /// Test timeConstant()
  void testTimeConstant(void);

  /// Test openTimeStep() and closeTimeStep()
  void testTimeStep(void);

  /// Test writeVertexField.
  void testWriteVertexField(void);

  /// Test writeCellField.
  void testWriteCellField(void);

  /// Test vtkFilename.
  void testVtkFilename(void);

  // PROTECTED MEMBERS //////////////////////////////////////////////////
protected :

  /// Initialize mesh.
  void _initialize(void);

  // PROTECTED MEMBERS //////////////////////////////////////////////////
protected :

  topology::Mesh* _mesh; ///< Mesh for data
  bool _flipFault; ///< If true, flip fault orientation.

  // PRIVATE MEMBERS ////////////////////////////////////////////////////
private :

  /** Create vertex fields.
   *
   * @param fields Vertex fields.
   */
  void
  _createVertexFields(topology::Fields<topology::Field<topology::Mesh> >* fields) const;

  /** Create cell fields.
   *
   * @param fields Cell fields.
   */
  void
  _createCellFields(topology::Fields<topology::Field<topology::Mesh> >* fields) const;

}; // class TestDataWriterVTKMesh

#endif // pylith_meshio_testdatawritervtkmesh_hh


// End of file 
