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
 * @file unittests/libtests/meshio/TestDataWriterVTKFaultMeshCases.hh
 *
 * @brief C++ unit testing for DataWriterVTK with fault mesh and
 * various cell types.
 */

#if !defined(pylith_meshio_testdatawritervtkfaultmeshcases_hh)
#define pylith_meshio_testdatawritervtkfaultmeshcases_hh

#include "TestDataWriterVTKFaultMesh.hh"

/// Namespace for pylith package
namespace pylith {
  namespace meshio {
    class TestDataWriterVTKFaultMeshTri3;
    class TestDataWriterVTKFaultMeshQuad4;
    class TestDataWriterVTKFaultMeshTet4;
    class TestDataWriterVTKFaultMeshHex8;
  } // meshio
} // pylith


// ----------------------------------------------------------------------
/// C++ unit testing for DataWriterVTK
class pylith::meshio::TestDataWriterVTKFaultMeshTri3 : public TestDataWriterVTKFaultMesh
{ // class TestDataWriterVTKFaultMeshTri3

  // CPPUNIT TEST SUITE /////////////////////////////////////////////////
  CPPUNIT_TEST_SUITE( TestDataWriterVTKFaultMeshTri3 );

  CPPUNIT_TEST( testTimeStep );
  CPPUNIT_TEST( testWriteVertexField );
  CPPUNIT_TEST( testWriteCellField );

  CPPUNIT_TEST_SUITE_END();

  // PUBLIC METHODS /////////////////////////////////////////////////////
public :

  /// Setup testing data.
  void setUp(void);

}; // class TestDataWriterVTKFaultMeshTri3


// ----------------------------------------------------------------------
/// C++ unit testing for DataWriterVTK
class pylith::meshio::TestDataWriterVTKFaultMeshQuad4 : public TestDataWriterVTKFaultMesh
{ // class TestDataWriterVTKFaultMeshQuad4

  // CPPUNIT TEST SUITE /////////////////////////////////////////////////
  CPPUNIT_TEST_SUITE( TestDataWriterVTKFaultMeshQuad4 );

  CPPUNIT_TEST( testTimeStep );
  CPPUNIT_TEST( testWriteVertexField );
  CPPUNIT_TEST( testWriteCellField );

  CPPUNIT_TEST_SUITE_END();

  // PUBLIC METHODS /////////////////////////////////////////////////////
public :

  /// Setup testing data.
  void setUp(void);

}; // class TestDataWriterVTKFaultMeshQuad4


// ----------------------------------------------------------------------
/// C++ unit testing for DataWriterVTK
class pylith::meshio::TestDataWriterVTKFaultMeshTet4 : public TestDataWriterVTKFaultMesh
{ // class TestDataWriterVTKFaultMeshTet4

  // CPPUNIT TEST SUITE /////////////////////////////////////////////////
  CPPUNIT_TEST_SUITE( TestDataWriterVTKFaultMeshTet4 );

  CPPUNIT_TEST( testTimeStep );
  CPPUNIT_TEST( testWriteVertexField );
  CPPUNIT_TEST( testWriteCellField );

  CPPUNIT_TEST_SUITE_END();

  // PUBLIC METHODS /////////////////////////////////////////////////////
public :

  /// Setup testing data.
  void setUp(void);

}; // class TestDataWriterVTKFaultMeshTet4


// ----------------------------------------------------------------------
/// C++ unit testing for DataWriterVTK
class pylith::meshio::TestDataWriterVTKFaultMeshHex8 : public TestDataWriterVTKFaultMesh
{ // class TestDataWriterVTKFaultMeshHex8

  // CPPUNIT TEST SUITE /////////////////////////////////////////////////
  CPPUNIT_TEST_SUITE( TestDataWriterVTKFaultMeshHex8 );

  CPPUNIT_TEST( testTimeStep );
  CPPUNIT_TEST( testWriteVertexField );
  CPPUNIT_TEST( testWriteCellField );

  CPPUNIT_TEST_SUITE_END();

  // PUBLIC METHODS /////////////////////////////////////////////////////
public :

  /// Setup testing data.
  void setUp(void);

}; // class TestDataWriterVTKFaultMeshHex8


#endif // pylith_meshio_testdatawritervtkfaultmeshcases_hh


// End of file 
