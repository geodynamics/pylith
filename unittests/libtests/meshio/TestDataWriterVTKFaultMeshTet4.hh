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
 * @file unittests/libtests/meshio/TestDataWriterVTKFaultMeshTet4.hh
 *
 * @brief C++ TestDataWriterVTKFaultMeshTet4 object
 *
 * C++ unit testing for DataWriterVTKFaultMeshTet4.
 */

#if !defined(pylith_meshio_testdatawritervtkfaultmeshtet4_hh)
#define pylith_meshio_testdatawritervtkfaultmeshtet4_hh

#include "TestDataWriterVTKFaultMesh.hh"

/// Namespace for pylith package
namespace pylith {
  namespace meshio {
    class TestDataWriterVTKFaultMeshTet4;
  } // meshio
} // pylith

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

#endif // pylith_meshio_testdatawritervtkfaultmeshtet4_hh


// End of file 
