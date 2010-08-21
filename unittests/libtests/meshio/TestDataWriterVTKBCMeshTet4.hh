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
// Copyright (c) 2010 University of California, Davis
//
// See COPYING for license information.
//
// ----------------------------------------------------------------------
//

/**
 * @file unittests/libtests/meshio/TestDataWriterVTKBCMeshTet4.hh
 *
 * @brief C++ TestDataWriterVTKBCMeshTet4 object
 *
 * C++ unit testing for DataWriterVTKBCMeshTet4.
 */

#if !defined(pylith_meshio_testdatawritervtkbcmeshtet4_hh)
#define pylith_meshio_testdatawritervtkbcmeshtet4_hh

#include "TestDataWriterVTKBCMesh.hh"

/// Namespace for pylith package
namespace pylith {
  namespace meshio {
    class TestDataWriterVTKBCMeshTet4;
  } // meshio
} // pylith

/// C++ unit testing for DataWriterVTK
class pylith::meshio::TestDataWriterVTKBCMeshTet4 : public TestDataWriterVTKBCMesh
{ // class TestDataWriterVTKBCMeshTet4

  // CPPUNIT TEST SUITE /////////////////////////////////////////////////
  CPPUNIT_TEST_SUITE( TestDataWriterVTKBCMeshTet4 );

  CPPUNIT_TEST( testTimeStep );
  CPPUNIT_TEST( testWriteVertexField );
  CPPUNIT_TEST( testWriteCellField );

  CPPUNIT_TEST_SUITE_END();

  // PUBLIC METHODS /////////////////////////////////////////////////////
public :

  /// Setup testing data.
  void setUp(void);

}; // class TestDataWriterVTKBCMeshTet4

#endif // pylith_meshio_testdatawritervtkbcmeshtet4_hh


// End of file 
