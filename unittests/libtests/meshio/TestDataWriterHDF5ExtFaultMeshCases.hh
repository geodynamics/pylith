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
 * @file unittests/libtests/meshio/TestDataWriterHDF5ExtFaultMeshCases.hh
 *
 * @brief C++ unit testing for DataWriterHDF5 with fault mesh and
 * various cell types.
 */

#if !defined(pylith_meshio_testdatawriterhdf5extfaultmeshcases_hh)
#define pylith_meshio_testdatawriterhdf5extfaultmeshcases_hh

#include "TestDataWriterHDF5ExtFaultMesh.hh"

/// Namespace for pylith package
namespace pylith {
  namespace meshio {
    class TestDataWriterHDF5ExtFaultMeshTri3;
    class TestDataWriterHDF5ExtFaultMeshQuad4;
    class TestDataWriterHDF5ExtFaultMeshTet4;
    class TestDataWriterHDF5ExtFaultMeshHex8;
  } // meshio
} // pylith


// ----------------------------------------------------------------------
/// C++ unit testing for DataWriterHDF5
class pylith::meshio::TestDataWriterHDF5ExtFaultMeshTri3 : public TestDataWriterHDF5ExtFaultMesh
{ // class TestDataWriterHDF5ExtFaultMeshTri3

  // CPPUNIT TEST SUITE /////////////////////////////////////////////////
  CPPUNIT_TEST_SUITE( TestDataWriterHDF5ExtFaultMeshTri3 );

  CPPUNIT_TEST( testOpenClose );
  CPPUNIT_TEST( testWriteVertexField );
  CPPUNIT_TEST( testWriteCellField );

  CPPUNIT_TEST_SUITE_END();

  // PUBLIC METHODS /////////////////////////////////////////////////////
public :

  /// Setup testing data.
  void setUp(void);

}; // class TestDataWriterHDF5ExtFaultMeshTri3


// ----------------------------------------------------------------------
/// C++ unit testing for DataWriterHDF5
class pylith::meshio::TestDataWriterHDF5ExtFaultMeshQuad4 : public TestDataWriterHDF5ExtFaultMesh
{ // class TestDataWriterHDF5ExtFaultMeshQuad4

  // CPPUNIT TEST SUITE /////////////////////////////////////////////////
  CPPUNIT_TEST_SUITE( TestDataWriterHDF5ExtFaultMeshQuad4 );

  CPPUNIT_TEST( testOpenClose );
  CPPUNIT_TEST( testWriteVertexField );
  CPPUNIT_TEST( testWriteCellField );

  CPPUNIT_TEST_SUITE_END();

  // PUBLIC METHODS /////////////////////////////////////////////////////
public :

  /// Setup testing data.
  void setUp(void);

}; // class TestDataWriterHDF5ExtFaultMeshQuad4


// ----------------------------------------------------------------------
/// C++ unit testing for DataWriterHDF5
class pylith::meshio::TestDataWriterHDF5ExtFaultMeshTet4 : public TestDataWriterHDF5ExtFaultMesh
{ // class TestDataWriterHDF5ExtFaultMeshTet4

  // CPPUNIT TEST SUITE /////////////////////////////////////////////////
  CPPUNIT_TEST_SUITE( TestDataWriterHDF5ExtFaultMeshTet4 );

  CPPUNIT_TEST( testOpenClose );
  CPPUNIT_TEST( testWriteVertexField );
  CPPUNIT_TEST( testWriteCellField );

  CPPUNIT_TEST_SUITE_END();

  // PUBLIC METHODS /////////////////////////////////////////////////////
public :

  /// Setup testing data.
  void setUp(void);

}; // class TestDataWriterHDF5ExtFaultMeshTet4


// ----------------------------------------------------------------------
/// C++ unit testing for DataWriterHDF5
class pylith::meshio::TestDataWriterHDF5ExtFaultMeshHex8 : public TestDataWriterHDF5ExtFaultMesh
{ // class TestDataWriterHDF5ExtFaultMeshHex8

  // CPPUNIT TEST SUITE /////////////////////////////////////////////////
  CPPUNIT_TEST_SUITE( TestDataWriterHDF5ExtFaultMeshHex8 );

  CPPUNIT_TEST( testOpenClose );
  CPPUNIT_TEST( testWriteVertexField );
  CPPUNIT_TEST( testWriteCellField );

  CPPUNIT_TEST_SUITE_END();

  // PUBLIC METHODS /////////////////////////////////////////////////////
public :

  /// Setup testing data.
  void setUp(void);

}; // class TestDataWriterHDF5ExtFaultMeshHex8


#endif // pylith_meshio_testdatawriterhdf5extfaultmeshcases_hh


// End of file 
