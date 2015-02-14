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
 * @file unittests/libtests/meshio/TestDataWriterHDF5ExtBCMeshCases.hh
 *
 * @brief C++ unit testing for DataWriterHDF5 with bc mesh output for
 * various cell types.
 */

#if !defined(pylith_meshio_testdatawriterhdf5extbcmeshcases_hh)
#define pylith_meshio_testdatawriterhdf5extbcmeshcases_hh

#include "TestDataWriterHDF5ExtBCMesh.hh"

/// Namespace for pylith package
namespace pylith {
  namespace meshio {
    class TestDataWriterHDF5ExtBCMeshTri3;
    class TestDataWriterHDF5ExtBCMeshQuad4;
    class TestDataWriterHDF5ExtBCMeshTet4;
    class TestDataWriterHDF5ExtBCMeshHex8;
  } // meshio
} // pylith


// ----------------------------------------------------------------------
/// C++ unit testing for DataWriterHDF5
class pylith::meshio::TestDataWriterHDF5ExtBCMeshTri3 : public TestDataWriterHDF5ExtBCMesh
{ // class TestDataWriterHDF5ExtBCMeshTri3

  // CPPUNIT TEST SUITE /////////////////////////////////////////////////
  CPPUNIT_TEST_SUITE( TestDataWriterHDF5ExtBCMeshTri3 );

  CPPUNIT_TEST( testOpenClose );
  CPPUNIT_TEST( testWriteVertexField );
  CPPUNIT_TEST( testWriteCellField );

  CPPUNIT_TEST_SUITE_END();

  // PUBLIC METHODS /////////////////////////////////////////////////////
public :

  /// Setup testing data.
  void setUp(void);

}; // class TestDataWriterHDF5ExtBCMeshTri3


// ----------------------------------------------------------------------
/// C++ unit testing for DataWriterHDF5
class pylith::meshio::TestDataWriterHDF5ExtBCMeshQuad4 : public TestDataWriterHDF5ExtBCMesh
{ // class TestDataWriterHDF5ExtBCMeshQuad4

  // CPPUNIT TEST SUITE /////////////////////////////////////////////////
  CPPUNIT_TEST_SUITE( TestDataWriterHDF5ExtBCMeshQuad4 );

  CPPUNIT_TEST( testOpenClose );
  CPPUNIT_TEST( testWriteVertexField );
  CPPUNIT_TEST( testWriteCellField );

  CPPUNIT_TEST_SUITE_END();

  // PUBLIC METHODS /////////////////////////////////////////////////////
public :

  /// Setup testing data.
  void setUp(void);

}; // class TestDataWriterHDF5ExtBCMeshQuad4


// ----------------------------------------------------------------------
/// C++ unit testing for DataWriterHDF5
class pylith::meshio::TestDataWriterHDF5ExtBCMeshTet4 : public TestDataWriterHDF5ExtBCMesh
{ // class TestDataWriterHDF5ExtBCMeshTet4

  // CPPUNIT TEST SUITE /////////////////////////////////////////////////
  CPPUNIT_TEST_SUITE( TestDataWriterHDF5ExtBCMeshTet4 );

  CPPUNIT_TEST( testOpenClose );
  CPPUNIT_TEST( testWriteVertexField );
  CPPUNIT_TEST( testWriteCellField );

  CPPUNIT_TEST_SUITE_END();

  // PUBLIC METHODS /////////////////////////////////////////////////////
public :

  /// Setup testing data.
  void setUp(void);

}; // class TestDataWriterHDF5ExtBCMeshTet4


// ----------------------------------------------------------------------
/// C++ unit testing for DataWriterHDF5
class pylith::meshio::TestDataWriterHDF5ExtBCMeshHex8 : public TestDataWriterHDF5ExtBCMesh
{ // class TestDataWriterHDF5ExtBCMeshHex8

  // CPPUNIT TEST SUITE /////////////////////////////////////////////////
  CPPUNIT_TEST_SUITE( TestDataWriterHDF5ExtBCMeshHex8 );

  CPPUNIT_TEST( testOpenClose );
  CPPUNIT_TEST( testWriteVertexField );
  CPPUNIT_TEST( testWriteCellField );

  CPPUNIT_TEST_SUITE_END();

  // PUBLIC METHODS /////////////////////////////////////////////////////
public :

  /// Setup testing data.
  void setUp(void);

}; // class TestDataWriterHDF5ExtBCMeshHex8


#endif // pylith_meshio_testdatawriterhdf5extbcmeshcases_hh


// End of file 
