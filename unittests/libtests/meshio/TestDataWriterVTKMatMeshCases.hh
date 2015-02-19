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
 * @file unittests/libtests/meshio/TestDataWriterVTKMatMeshCases.hh
 *
 * @brief C++ unit testing for DataWriterVTK with material mesh and
 * various cell types.
 */

#if !defined(pylith_meshio_testdatawritervtkmatmeshcases_hh)
#define pylith_meshio_testdatawritervtkmatmeshcases_hh

#include "TestDataWriterVTKMesh.hh"

/// Namespace for pylith package
namespace pylith {
  namespace meshio {
    class TestDataWriterVTKMatMeshTri3;
    class TestDataWriterVTKMatMeshQuad4;
    class TestDataWriterVTKMatMeshTet4;
    class TestDataWriterVTKMatMeshHex8;
  } // meshio
} // pylith


// ----------------------------------------------------------------------
/// C++ unit testing for DataWriterVTK
class pylith::meshio::TestDataWriterVTKMatMeshTri3 : public TestDataWriterVTKMesh
{ // class TestDataWriterVTKMatMeshTri3

  // CPPUNIT TEST SUITE /////////////////////////////////////////////////
  CPPUNIT_TEST_SUITE( TestDataWriterVTKMatMeshTri3 );

  CPPUNIT_TEST( testTimeStep );
  CPPUNIT_TEST( testWriteVertexField );
  CPPUNIT_TEST( testWriteCellField );

  CPPUNIT_TEST_SUITE_END();

  // PUBLIC METHODS /////////////////////////////////////////////////////
public :

  /// Setup testing data.
  void setUp(void);

}; // class TestDataWriterVTKMatMeshTri3


// ----------------------------------------------------------------------
/// C++ unit testing for DataWriterVTK
class pylith::meshio::TestDataWriterVTKMatMeshQuad4 : public TestDataWriterVTKMesh
{ // class TestDataWriterVTKMatMeshQuad4

  // CPPUNIT TEST SUITE /////////////////////////////////////////////////
  CPPUNIT_TEST_SUITE( TestDataWriterVTKMatMeshQuad4 );

  CPPUNIT_TEST( testTimeStep );
  CPPUNIT_TEST( testWriteVertexField );
  CPPUNIT_TEST( testWriteCellField );

  CPPUNIT_TEST_SUITE_END();

  // PUBLIC METHODS /////////////////////////////////////////////////////
public :

  /// Setup testing data.
  void setUp(void);

}; // class TestDataWriterVTKMatMeshQuad4


// ----------------------------------------------------------------------
/// C++ unit testing for DataWriterVTK
class pylith::meshio::TestDataWriterVTKMatMeshTet4 : public TestDataWriterVTKMesh
{ // class TestDataWriterVTKMatMeshTet4

  // CPPUNIT TEST SUITE /////////////////////////////////////////////////
  CPPUNIT_TEST_SUITE( TestDataWriterVTKMatMeshTet4 );

  CPPUNIT_TEST( testTimeStep );
  CPPUNIT_TEST( testWriteVertexField );
  CPPUNIT_TEST( testWriteCellField );

  CPPUNIT_TEST_SUITE_END();

  // PUBLIC METHODS /////////////////////////////////////////////////////
public :

  /// Setup testing data.
  void setUp(void);

}; // class TestDataWriterVTKMatMeshTet4


// ----------------------------------------------------------------------
/// C++ unit testing for DataWriterVTK
class pylith::meshio::TestDataWriterVTKMatMeshHex8 : public TestDataWriterVTKMesh
{ // class TestDataWriterVTKMatMeshHex8

  // CPPUNIT TEST SUITE /////////////////////////////////////////////////
  CPPUNIT_TEST_SUITE( TestDataWriterVTKMatMeshHex8 );

  CPPUNIT_TEST( testTimeStep );
  CPPUNIT_TEST( testWriteVertexField );
  CPPUNIT_TEST( testWriteCellField );

  CPPUNIT_TEST_SUITE_END();

  // PUBLIC METHODS /////////////////////////////////////////////////////
public :

  /// Setup testing data.
  void setUp(void);

}; // class TestDataWriterVTKMatMeshHex8


#endif // pylith_meshio_testdatawritervtkmatmeshcases_hh


// End of file 
