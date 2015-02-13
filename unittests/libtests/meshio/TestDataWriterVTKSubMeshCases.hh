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
 * @file unittests/libtests/meshio/TestDataWriterVTKSubMeshCases.hh
 *
 * @brief C++ unit testing for DataWriterVTK with subdomain mesh and
 * various cell types.
 */

#if !defined(pylith_meshio_testdatawritervtksubmeshcases_hh)
#define pylith_meshio_testdatawritervtksubmeshcases_hh

#include "TestDataWriterVTKSubMesh.hh"

/// Namespace for pylith package
namespace pylith {
  namespace meshio {
    class TestDataWriterVTKSubMeshTri3;
    class TestDataWriterVTKSubMeshQuad4;
    class TestDataWriterVTKSubMeshTet4;
    class TestDataWriterVTKSubMeshHex8;
  } // meshio
} // pylith


// ----------------------------------------------------------------------
/// C++ unit testing for DataWriterVTK
class pylith::meshio::TestDataWriterVTKSubMeshTri3 : public TestDataWriterVTKSubMesh
{ // class TestDataWriterVTKSubMeshTri3

  // CPPUNIT TEST SUITE /////////////////////////////////////////////////
  CPPUNIT_TEST_SUITE( TestDataWriterVTKSubMeshTri3 );

  CPPUNIT_TEST( testTimeStep );
  CPPUNIT_TEST( testWriteVertexField );
  CPPUNIT_TEST( testWriteCellField );

  CPPUNIT_TEST_SUITE_END();

  // PUBLIC METHODS /////////////////////////////////////////////////////
public :

  /// Setup testing data.
  void setUp(void);

}; // class TestDataWriterVTKSubMeshTri3


// ----------------------------------------------------------------------
/// C++ unit testing for DataWriterVTK
class pylith::meshio::TestDataWriterVTKSubMeshQuad4 : public TestDataWriterVTKSubMesh
{ // class TestDataWriterVTKSubMeshQuad4

  // CPPUNIT TEST SUITE /////////////////////////////////////////////////
  CPPUNIT_TEST_SUITE( TestDataWriterVTKSubMeshQuad4 );

  CPPUNIT_TEST( testTimeStep );
  CPPUNIT_TEST( testWriteVertexField );
  CPPUNIT_TEST( testWriteCellField );

  CPPUNIT_TEST_SUITE_END();

  // PUBLIC METHODS /////////////////////////////////////////////////////
public :

  /// Setup testing data.
  void setUp(void);

}; // class TestDataWriterVTKSubMeshQuad4


// ----------------------------------------------------------------------
/// C++ unit testing for DataWriterVTK
class pylith::meshio::TestDataWriterVTKSubMeshTet4 : public TestDataWriterVTKSubMesh
{ // class TestDataWriterVTKSubMeshTet4

  // CPPUNIT TEST SUITE /////////////////////////////////////////////////
  CPPUNIT_TEST_SUITE( TestDataWriterVTKSubMeshTet4 );

  CPPUNIT_TEST( testTimeStep );
  CPPUNIT_TEST( testWriteVertexField );
  CPPUNIT_TEST( testWriteCellField );

  CPPUNIT_TEST_SUITE_END();

  // PUBLIC METHODS /////////////////////////////////////////////////////
public :

  /// Setup testing data.
  void setUp(void);

}; // class TestDataWriterVTKSubMeshTet4


// ----------------------------------------------------------------------
/// C++ unit testing for DataWriterVTK
class pylith::meshio::TestDataWriterVTKSubMeshHex8 : public TestDataWriterVTKSubMesh
{ // class TestDataWriterVTKSubMeshHex8

  // CPPUNIT TEST SUITE /////////////////////////////////////////////////
  CPPUNIT_TEST_SUITE( TestDataWriterVTKSubMeshHex8 );

  CPPUNIT_TEST( testTimeStep );
  CPPUNIT_TEST( testWriteVertexField );
  CPPUNIT_TEST( testWriteCellField );

  CPPUNIT_TEST_SUITE_END();

  // PUBLIC METHODS /////////////////////////////////////////////////////
public :

  /// Setup testing data.
  void setUp(void);

}; // class TestDataWriterVTKSubMeshHex8


#endif // pylith_meshio_testdatawritervtksubmeshcases_hh


// End of file 
