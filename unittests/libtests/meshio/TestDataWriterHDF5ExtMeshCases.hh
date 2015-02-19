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
 * @file unittests/libtests/meshio/TestDataWriterHDF5ExtMeshCases.hh
 *
 * @brief C++ unit testing for DataWriterHDF5Ext with mesh and various
 * cell types.
 */

#if !defined(pylith_meshio_testdatawriterhdf5extmeshcases_hh)
#define pylith_meshio_testdatawriterhdf5extmeshcases_hh

#include "TestDataWriterHDF5ExtMesh.hh"

/// Namespace for pylith package
namespace pylith {
  namespace meshio {
    class TestDataWriterHDF5ExtMeshTri3;
    class TestDataWriterHDF5ExtMeshQuad4;
    class TestDataWriterHDF5ExtMeshTet4;
    class TestDataWriterHDF5ExtMeshHex8;
  } // meshio
} // pylith


// ----------------------------------------------------------------------
/// C++ unit testing for DataWriterHDF5
class pylith::meshio::TestDataWriterHDF5ExtMeshTri3 : public TestDataWriterHDF5ExtMesh
{ // class TestDataWriterHDF5ExtMeshTri3

  // CPPUNIT TEST SUITE /////////////////////////////////////////////////
  CPPUNIT_TEST_SUITE( TestDataWriterHDF5ExtMeshTri3 );

  CPPUNIT_TEST( testOpenClose );
  CPPUNIT_TEST( testWriteVertexField );
  CPPUNIT_TEST( testWriteCellField );

  CPPUNIT_TEST_SUITE_END();

  // PUBLIC METHODS /////////////////////////////////////////////////////
public :

  /// Setup testing data.
  void setUp(void);

}; // class TestDataWriterHDF5ExtMeshTri3


// ----------------------------------------------------------------------
/// C++ unit testing for DataWriterHDF5
class pylith::meshio::TestDataWriterHDF5ExtMeshQuad4 : public TestDataWriterHDF5ExtMesh
{ // class TestDataWriterHDF5ExtMeshQuad4

  // CPPUNIT TEST SUITE /////////////////////////////////////////////////
  CPPUNIT_TEST_SUITE( TestDataWriterHDF5ExtMeshQuad4 );

  CPPUNIT_TEST( testOpenClose );
  CPPUNIT_TEST( testWriteVertexField );
  CPPUNIT_TEST( testWriteCellField );

  CPPUNIT_TEST_SUITE_END();

  // PUBLIC METHODS /////////////////////////////////////////////////////
public :

  /// Setup testing data.
  void setUp(void);

}; // class TestDataWriterHDF5ExtMeshQuad4


// ----------------------------------------------------------------------
/// C++ unit testing for DataWriterHDF5
class pylith::meshio::TestDataWriterHDF5ExtMeshTet4 : public TestDataWriterHDF5ExtMesh
{ // class TestDataWriterHDF5ExtMeshTet4

  // CPPUNIT TEST SUITE /////////////////////////////////////////////////
  CPPUNIT_TEST_SUITE( TestDataWriterHDF5ExtMeshTet4 );

  CPPUNIT_TEST( testOpenClose );
  CPPUNIT_TEST( testWriteVertexField );
  CPPUNIT_TEST( testWriteCellField );

  CPPUNIT_TEST_SUITE_END();

  // PUBLIC METHODS /////////////////////////////////////////////////////
public :

  /// Setup testing data.
  void setUp(void);

}; // class TestDataWriterHDF5ExtMeshTet4


// ----------------------------------------------------------------------
/// C++ unit testing for DataWriterHDF5
class pylith::meshio::TestDataWriterHDF5ExtMeshHex8 : public TestDataWriterHDF5ExtMesh
{ // class TestDataWriterHDF5ExtMeshHex8

  // CPPUNIT TEST SUITE /////////////////////////////////////////////////
  CPPUNIT_TEST_SUITE( TestDataWriterHDF5ExtMeshHex8 );

  CPPUNIT_TEST( testOpenClose );
  CPPUNIT_TEST( testWriteVertexField );
  CPPUNIT_TEST( testWriteCellField );

  CPPUNIT_TEST_SUITE_END();

  // PUBLIC METHODS /////////////////////////////////////////////////////
public :

  /// Setup testing data.
  void setUp(void);

}; // class TestDataWriterHDF5ExtMeshHex8


#endif // pylith_meshio_testdatawriterhdf5extmeshcases_hh


// End of file 
