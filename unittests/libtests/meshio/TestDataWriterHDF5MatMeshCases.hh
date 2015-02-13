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
 * @file unittests/libtests/meshio/TestDataWriterHDF5MatMeshCases.hh
 *
 * @brief C++ unit testing for DataWriterHDF5 with material mesh and
 * various cell types.
 */

#if !defined(pylith_meshio_testdatawriterhdf5matmeshcases_hh)
#define pylith_meshio_testdatawriterhdf5matmeshcases_hh

#include "TestDataWriterHDF5Mesh.hh"

/// Namespace for pylith package
namespace pylith {
  namespace meshio {
    class TestDataWriterHDF5MatMeshTri3;
    class TestDataWriterHDF5MatMeshQuad4;
    class TestDataWriterHDF5MatMeshTet4;
    class TestDataWriterHDF5MatMeshHex8;
  } // meshio
} // pylith


// ----------------------------------------------------------------------
/// C++ unit testing for DataWriterHDF5
class pylith::meshio::TestDataWriterHDF5MatMeshTri3 : public TestDataWriterHDF5Mesh
{ // class TestDataWriterHDF5MatMeshTri3

  // CPPUNIT TEST SUITE /////////////////////////////////////////////////
  CPPUNIT_TEST_SUITE( TestDataWriterHDF5MatMeshTri3 );

  CPPUNIT_TEST( testOpenClose );
  CPPUNIT_TEST( testWriteVertexField );
  CPPUNIT_TEST( testWriteCellField );

  CPPUNIT_TEST_SUITE_END();

  // PUBLIC METHODS /////////////////////////////////////////////////////
public :

  /// Setup testing data.
  void setUp(void);

}; // class TestDataWriterHDF5MatMeshTri3


// ----------------------------------------------------------------------
/// C++ unit testing for DataWriterHDF5
class pylith::meshio::TestDataWriterHDF5MatMeshQuad4 : public TestDataWriterHDF5Mesh
{ // class TestDataWriterHDF5MatMeshQuad4

  // CPPUNIT TEST SUITE /////////////////////////////////////////////////
  CPPUNIT_TEST_SUITE( TestDataWriterHDF5MatMeshQuad4 );

  CPPUNIT_TEST( testOpenClose );
  CPPUNIT_TEST( testWriteVertexField );
  CPPUNIT_TEST( testWriteCellField );

  CPPUNIT_TEST_SUITE_END();

  // PUBLIC METHODS /////////////////////////////////////////////////////
public :

  /// Setup testing data.
  void setUp(void);

}; // class TestDataWriterHDF5MatMeshQuad4


// ----------------------------------------------------------------------
/// C++ unit testing for DataWriterHDF5
class pylith::meshio::TestDataWriterHDF5MatMeshTet4 : public TestDataWriterHDF5Mesh
{ // class TestDataWriterHDF5MatMeshTet4

  // CPPUNIT TEST SUITE /////////////////////////////////////////////////
  CPPUNIT_TEST_SUITE( TestDataWriterHDF5MatMeshTet4 );

  CPPUNIT_TEST( testOpenClose );
  CPPUNIT_TEST( testWriteVertexField );
  CPPUNIT_TEST( testWriteCellField );

  CPPUNIT_TEST_SUITE_END();

  // PUBLIC METHODS /////////////////////////////////////////////////////
public :

  /// Setup testing data.
  void setUp(void);

}; // class TestDataWriterHDF5MatMeshTet4


// ----------------------------------------------------------------------
/// C++ unit testing for DataWriterHDF5
class pylith::meshio::TestDataWriterHDF5MatMeshHex8 : public TestDataWriterHDF5Mesh
{ // class TestDataWriterHDF5MatMeshHex8

  // CPPUNIT TEST SUITE /////////////////////////////////////////////////
  CPPUNIT_TEST_SUITE( TestDataWriterHDF5MatMeshHex8 );

  CPPUNIT_TEST( testOpenClose );
  CPPUNIT_TEST( testWriteVertexField );
  CPPUNIT_TEST( testWriteCellField );

  CPPUNIT_TEST_SUITE_END();

  // PUBLIC METHODS /////////////////////////////////////////////////////
public :

  /// Setup testing data.
  void setUp(void);

}; // class TestDataWriterHDF5MatMeshHex8


#endif // pylith_meshio_testdatawriterhdf5matmeshcases_hh


// End of file 
