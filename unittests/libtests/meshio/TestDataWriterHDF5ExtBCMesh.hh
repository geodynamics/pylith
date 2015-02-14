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
 * @file unittests/libtests/meshio/TestDataWriterHDF5ExtBCMesh.hh
 *
 * @brief C++ TestDataWriterHDF5ExtBCMesh object
 *
 * C++ unit testing for DataWriterHDF5ExtBCMesh.
 */

#if !defined(pylith_meshio_testdatawriterhdf5extbcmesh_hh)
#define pylith_meshio_testdatawriterhdf5extbcmesh_hh

#include "TestDataWriterHDF5.hh"
#include "TestDataWriterBCMesh.hh"

#include "pylith/topology/topologyfwd.hh" // USES Mesh, Field

#include <cppunit/extensions/HelperMacros.h>

/// Namespace for pylith package
namespace pylith {
  namespace meshio {
    class TestDataWriterHDF5ExtBCMesh;
  } // meshio
} // pylith

/// C++ unit testing for DataWriterHDF5
class pylith::meshio::TestDataWriterHDF5ExtBCMesh : public TestDataWriterHDF5,
						    public TestDataWriterBCMesh,
						    public CppUnit::TestFixture
{ // class TestDataWriterHDF5ExtBCMesh

  // CPPUNIT TEST SUITE /////////////////////////////////////////////////
  CPPUNIT_TEST_SUITE( TestDataWriterHDF5ExtBCMesh );

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

  /// Test open() and close().
  void testOpenClose(void);

  /// Test writeVertexField.
  void testWriteVertexField(void);

  /// Test writeCellField.
  void testWriteCellField(void);

}; // class TestDataWriterHDF5ExtBCMesh

#endif // pylith_meshio_testdatawriterhdf5extbcmesh_hh


// End of file 
