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
 * @file unittests/libtests/meshio/TestDataWriterHDF5ExtSubMesh.hh
 *
 * @brief C++ TestDataWriterHDF5ExtSubMesh object
 *
 * C++ unit testing for DataWriterHDF5ExtSubMesh.
 */

#if !defined(pylith_meshio_testdatawriterhdf5extsubmesh_hh)
#define pylith_meshio_testdatawriterhdf5extsubmesh_hh

#include "TestDataWriterHDF5.hh"
#include "TestDataWriterSubMesh.hh"

#include "pylith/topology/topologyfwd.hh" // USES Mesh, Field

#include <cppunit/extensions/HelperMacros.h>

/// Namespace for pylith package
namespace pylith {
  namespace meshio {
    class TestDataWriterHDF5ExtSubMesh;
  } // meshio
} // pylith

/// C++ unit testing for DataWriterHDF5
class pylith::meshio::TestDataWriterHDF5ExtSubMesh : public TestDataWriterHDF5,
						     public TestDataWriterSubMesh,
						     public CppUnit::TestFixture
{ // class TestDataWriterHDF5ExtSubMesh

  // CPPUNIT TEST SUITE /////////////////////////////////////////////////
  CPPUNIT_TEST_SUITE( TestDataWriterHDF5ExtSubMesh );

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

  /// Test open() and close()
  void testOpenClose(void);

  /// Test writeVertexField.
  void testWriteVertexField(void);

  /// Test writeCellField.
  void testWriteCellField(void);

}; // class TestDataWriterHDF5ExtSubMesh

#endif // pylith_meshio_testdatawriterhdf5extsubmesh_hh


// End of file 
