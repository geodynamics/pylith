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
 * @file unittests/libtests/meshio/TestDataWriterVTK.hh
 *
 * @brief C++ TestDataWriterVTK object
 *
 * C++ unit testing for DataWriterVTK.
 */

#if !defined(pylith_meshio_testdatawritervtk_hh)
#define pylith_meshio_testdatawritervtk_hh

#include "pylith/topology/topologyfwd.hh" // USES Mesh, Field

#include <cppunit/extensions/HelperMacros.h>

/// Namespace for pylith package
namespace pylith {
  namespace meshio {
    class TestDataWriterVTK;

    class DataWriterVTKData;
  } // meshio
} // pylith

/// C++ unit testing for DataWriterVTK
class pylith::meshio::TestDataWriterVTK : public CppUnit::TestFixture
{ // class TestDataWriterVTK

  // CPPUNIT TEST SUITE /////////////////////////////////////////////////
  CPPUNIT_TEST_SUITE( TestDataWriterVTK );

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

  // PUBLIC METHODS /////////////////////////////////////////////////////
public :

  /** Check VTK file against archived file.
   *
   * @param filename Name of file to check.
   * @param t Time for file.
   * @param timeFormat Format of timestamp in filename.
   */
  static
  void checkFile(const char* filename,
		 const double t,
		 const char* timeFormat);
  
  // PROTECTED MEMBERS //////////////////////////////////////////////////
protected :

  DataWriterVTKData* _data; ///< Data for testing
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

}; // class TestDataWriterVTK

#endif // pylith_meshio_testdatawritervtk_hh


// End of file 
