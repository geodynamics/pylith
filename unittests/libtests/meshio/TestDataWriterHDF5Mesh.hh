// -*- C++ -*-
//
// ----------------------------------------------------------------------
//
// Brad T. Aagaard, U.S. Geological Survey
// Charles A. Williams, GNS Science
// Matthew G. Knepley, University of Chicago
//
// This code was developed as part of the Computational Infrastructure
// for Geodynamics (http://geodynamics.org).
//
// Copyright (c) 2010-2015 University of California, Davis
//
// See COPYING for license information.
//
// ----------------------------------------------------------------------
//

/**
 * @file unittests/libtests/meshio/TestDataWriterHDF5Mesh.hh
 *
 * @brief C++ TestDataWriterHDF5Mesh object
 *
 * C++ unit testing for DataWriterHDF5Mesh.
 */

#if !defined(pylith_meshio_testdatawriterhdf5mesh_hh)
#define pylith_meshio_testdatawriterhdf5mesh_hh

#include "TestDataWriterHDF5.hh" // ISA TestDataWriterHDF5
#include "TestDataWriterMesh.hh" // ISA TestDataWriterMesh

#include <cppunit/extensions/HelperMacros.h>

/// Namespace for pylith package
namespace pylith {
  namespace meshio {
    class TestDataWriterHDF5Mesh;
  } // meshio
} // pylith

/// C++ unit testing for DataWriterHDF5
class pylith::meshio::TestDataWriterHDF5Mesh : public TestDataWriterHDF5,
					       public TestDataWriterMesh,
					       public CppUnit::TestFixture
{ // class TestDataWriterHDF5Mesh

  // CPPUNIT TEST SUITE /////////////////////////////////////////////////
  CPPUNIT_TEST_SUITE( TestDataWriterHDF5Mesh );

  CPPUNIT_TEST( testConstructor );
  CPPUNIT_TEST( testFilename );
  CPPUNIT_TEST( testHdf5Filename );

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

  /// Test open() and close()
  void testOpenClose(void);

  /// Test writeVertexField.
  void testWriteVertexField(void);

  /// Test writeCellField.
  void testWriteCellField(void);

  /// Test hdf5Filename.
  void testHdf5Filename(void);

}; // class TestDataWriterHDF5Mesh

#endif // pylith_meshio_testdatawriterhdf5mesh_hh


// End of file 
