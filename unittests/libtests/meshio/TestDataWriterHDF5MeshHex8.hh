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
 * @file unittests/libtests/meshio/TestDataWriterHDF5MeshHex8.hh
 *
 * @brief C++ TestDataWriterHDF5MeshHex8 object
 *
 * C++ unit testing for DataWriterHDF5MeshHex8.
 */

#if !defined(pylith_meshio_testdatawriterhdf5meshhex8_hh)
#define pylith_meshio_testdatawriterhdf5meshhex8_hh

#include "TestDataWriterHDF5Mesh.hh"

/// Namespace for pylith package
namespace pylith {
  namespace meshio {
    class TestDataWriterHDF5MeshHex8;
  } // meshio
} // pylith

/// C++ unit testing for DataWriterHDF5
class pylith::meshio::TestDataWriterHDF5MeshHex8 : public TestDataWriterHDF5Mesh
{ // class TestDataWriterHDF5MeshHex8

  // CPPUNIT TEST SUITE /////////////////////////////////////////////////
  CPPUNIT_TEST_SUITE( TestDataWriterHDF5MeshHex8 );

  CPPUNIT_TEST( testTimeStep );
  CPPUNIT_TEST( testWriteVertexField );
  CPPUNIT_TEST( testWriteCellField );

  CPPUNIT_TEST_SUITE_END();

  // PUBLIC METHODS /////////////////////////////////////////////////////
public :

  /// Setup testing data.
  void setUp(void);

}; // class TestDataWriterHDF5MeshHex8

#endif // pylith_meshio_testdatawriterhdf5meshhex8_hh


// End of file 
