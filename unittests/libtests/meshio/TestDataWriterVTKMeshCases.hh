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
 * @file unittests/libtests/meshio/TestDataWriterVTKMeshCases.hh
 *
 * @brief C++ unit testing for DataWriterVTK mesh output with various
 * cell types.
 */

#if !defined(pylith_meshio_testdatawritervtkmeshcases_hh)
#define pylith_meshio_testdatawritervtkmeshcases_hh

#include "TestDataWriterVTKMesh.hh"

/// Namespace for pylith package
namespace pylith {
  namespace meshio {
    class TestDataWriterVTKMeshTri3;
    class TestDataWriterVTKMeshQuad4;
    class TestDataWriterVTKMeshTet4;
    class TestDataWriterVTKMeshHex8;
  } // meshio
} // pylith


// ----------------------------------------------------------------------
/// C++ unit testing for DataWriterVTK
class pylith::meshio::TestDataWriterVTKMeshTri3 : public TestDataWriterVTKMesh
{ // class TestDataWriterVTKMeshTri3

  // CPPUNIT TEST SUITE /////////////////////////////////////////////////
  CPPUNIT_TEST_SUITE( TestDataWriterVTKMeshTri3 );

  CPPUNIT_TEST( testTimeStep );
  CPPUNIT_TEST( testWriteVertexField );
  CPPUNIT_TEST( testWriteCellField );

  CPPUNIT_TEST_SUITE_END();

  // PUBLIC METHODS /////////////////////////////////////////////////////
public :

  /// Setup testing data.
  void setUp(void);

}; // class TestDataWriterVTKMeshTri3


// ----------------------------------------------------------------------
/// C++ unit testing for DataWriterVTK
class pylith::meshio::TestDataWriterVTKMeshQuad4 : public TestDataWriterVTKMesh
{ // class TestDataWriterVTKMeshQuad4

  // CPPUNIT TEST SUITE /////////////////////////////////////////////////
  CPPUNIT_TEST_SUITE( TestDataWriterVTKMeshQuad4 );

  CPPUNIT_TEST( testTimeStep );
  CPPUNIT_TEST( testWriteVertexField );
  CPPUNIT_TEST( testWriteCellField );

  CPPUNIT_TEST_SUITE_END();

  // PUBLIC METHODS /////////////////////////////////////////////////////
public :

  /// Setup testing data.
  void setUp(void);

}; // class TestDataWriterVTKMeshQuad4


// ----------------------------------------------------------------------
/// C++ unit testing for DataWriterVTK
class pylith::meshio::TestDataWriterVTKMeshTet4 : public TestDataWriterVTKMesh
{ // class TestDataWriterVTKMeshTet4

  // CPPUNIT TEST SUITE /////////////////////////////////////////////////
  CPPUNIT_TEST_SUITE( TestDataWriterVTKMeshTet4 );

  CPPUNIT_TEST( testTimeStep );
  CPPUNIT_TEST( testWriteVertexField );
  CPPUNIT_TEST( testWriteCellField );

  CPPUNIT_TEST_SUITE_END();

  // PUBLIC METHODS /////////////////////////////////////////////////////
public :

  /// Setup testing data.
  void setUp(void);

}; // class TestDataWriterVTKMeshTet4


// ----------------------------------------------------------------------
/// C++ unit testing for DataWriterVTK
class pylith::meshio::TestDataWriterVTKMeshHex8 : public TestDataWriterVTKMesh
{ // class TestDataWriterVTKMeshHex8

  // CPPUNIT TEST SUITE /////////////////////////////////////////////////
  CPPUNIT_TEST_SUITE( TestDataWriterVTKMeshHex8 );

  CPPUNIT_TEST( testTimeStep );
  CPPUNIT_TEST( testWriteVertexField );
  CPPUNIT_TEST( testWriteCellField );

  CPPUNIT_TEST_SUITE_END();

  // PUBLIC METHODS /////////////////////////////////////////////////////
public :

  /// Setup testing data.
  void setUp(void);

}; // class TestDataWriterVTKMeshHex8


#endif // pylith_meshio_testdatawritervtkmeshcases_hh


// End of file 
