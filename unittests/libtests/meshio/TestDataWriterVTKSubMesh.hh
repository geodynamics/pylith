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
 * @file unittests/libtests/meshio/TestDataWriterVTKSubMesh.hh
 *
 * @brief C++ TestDataWriterVTKSubMesh object
 *
 * C++ unit testing for DataWriterVTKSubMesh.
 */

#if !defined(pylith_meshio_testdatawritervtksubmesh_hh)
#define pylith_meshio_testdatawritervtksubmesh_hh

#include "TestDataWriterVTK.hh"
#include "TestDataWriterSubMesh.hh"

#include "pylith/topology/topologyfwd.hh" // USES Mesh, Field

#include <cppunit/extensions/HelperMacros.h>

/// Namespace for pylith package
namespace pylith {
  namespace meshio {
    class TestDataWriterVTKSubMesh;
  } // meshio
} // pylith

/// C++ unit testing for DataWriterVTK
class pylith::meshio::TestDataWriterVTKSubMesh : public TestDataWriterVTK,
						 public TestDataWriterSubMesh,
						 public CppUnit::TestFixture
{ // class TestDataWriterVTKSubMesh

  // CPPUNIT TEST SUITE /////////////////////////////////////////////////
  CPPUNIT_TEST_SUITE( TestDataWriterVTKSubMesh );

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

}; // class TestDataWriterVTKSubMesh

#endif // pylith_meshio_testdatawritervtksubmesh_hh


// End of file 
