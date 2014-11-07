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
// Copyright (c) 2010-2014 University of California, Davis
//
// See COPYING for license information.
//
// ----------------------------------------------------------------------
//

/**
 * @file unittests/libtests/meshio/TestDataWriterHDF5ExtPointsCases.hh
 *
 * @brief C++ unit testing for DataWriterHDF5Ext output interpolated to
 * points with various mesh cell types.
 */

#if !defined(pylith_meshio_testdatawriterhdf5extpointscases_hh)
#define pylith_meshio_testdatawriterhdf5extpointscases_hh

#include "TestDataWriterHDF5ExtPoints.hh"

/// Namespace for pylith package
namespace pylith {
  namespace meshio {
    class TestDataWriterHDF5ExtPointsTri3;
    class TestDataWriterHDF5ExtPointsQuad4;
    class TestDataWriterHDF5ExtPointsTet4;
    class TestDataWriterHDF5ExtPointsHex8;
  } // meshio
} // pylith


// ----------------------------------------------------------------------
/// C++ unit testing for DataWriterHDF5Ext
class pylith::meshio::TestDataWriterHDF5ExtPointsTri3 : public TestDataWriterHDF5ExtPoints
{ // class TestDataWriterHDF5ExtPointsTri3

  // CPPUNIT TEST SUITE /////////////////////////////////////////////////
  CPPUNIT_TEST_SUITE( TestDataWriterHDF5ExtPointsTri3 );

  CPPUNIT_TEST( testTimeStep );
  CPPUNIT_TEST( testWriteVertexField );

  CPPUNIT_TEST_SUITE_END();

  // PUBLIC METHODS /////////////////////////////////////////////////////
public :

  /// Setup testing data.
  void setUp(void);

}; // class TestDataWriterHDF5ExtPointsTri3


// ----------------------------------------------------------------------
/// C++ unit testing for DataWriterHDF5Ext
class pylith::meshio::TestDataWriterHDF5ExtPointsQuad4 : public TestDataWriterHDF5ExtPoints
{ // class TestDataWriterHDF5ExtPointsQuad4

  // CPPUNIT TEST SUITE /////////////////////////////////////////////////
  CPPUNIT_TEST_SUITE( TestDataWriterHDF5ExtPointsQuad4 );

  CPPUNIT_TEST( testTimeStep );
  CPPUNIT_TEST( testWriteVertexField );

  CPPUNIT_TEST_SUITE_END();

  // PUBLIC METHODS /////////////////////////////////////////////////////
public :

  /// Setup testing data.
  void setUp(void);

}; // class TestDataWriterHDF5ExtPointsQuad4


// ----------------------------------------------------------------------
/// C++ unit testing for DataWriterHDF5Ext
class pylith::meshio::TestDataWriterHDF5ExtPointsTet4 : public TestDataWriterHDF5ExtPoints
{ // class TestDataWriterHDF5ExtPointsTet4

  // CPPUNIT TEST SUITE /////////////////////////////////////////////////
  CPPUNIT_TEST_SUITE( TestDataWriterHDF5ExtPointsTet4 );

  CPPUNIT_TEST( testTimeStep );
  CPPUNIT_TEST( testWriteVertexField );

  CPPUNIT_TEST_SUITE_END();

  // PUBLIC METHODS /////////////////////////////////////////////////////
public :

  /// Setup testing data.
  void setUp(void);

}; // class TestDataWriterHDF5ExtPointsTet4


// ----------------------------------------------------------------------
/// C++ unit testing for DataWriterHDF5Ext
class pylith::meshio::TestDataWriterHDF5ExtPointsHex8 : public TestDataWriterHDF5ExtPoints
{ // class TestDataWriterHDF5ExtPointsHex8

  // CPPUNIT TEST SUITE /////////////////////////////////////////////////
  CPPUNIT_TEST_SUITE( TestDataWriterHDF5ExtPointsHex8 );

  CPPUNIT_TEST( testTimeStep );
  CPPUNIT_TEST( testWriteVertexField );

  CPPUNIT_TEST_SUITE_END();

  // PUBLIC METHODS /////////////////////////////////////////////////////
public :

  /// Setup testing data.
  void setUp(void);

}; // class TestDataWriterHDF5ExtPointsHex8


#endif // pylith_meshio_testdatawriterhdf5extpointscases_hh


// End of file 
