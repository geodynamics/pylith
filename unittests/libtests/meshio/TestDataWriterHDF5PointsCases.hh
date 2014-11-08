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
 * @file unittests/libtests/meshio/TestDataWriterHDF5PointsCases.hh
 *
 * @brief C++ unit testing for DataWriterHDF5 output interpolated to
 * points with various mesh cell types.
 */

#if !defined(pylith_meshio_testdatawriterhdf5pointscases_hh)
#define pylith_meshio_testdatawriterhdf5pointscases_hh

#include "TestDataWriterHDF5Points.hh"

/// Namespace for pylith package
namespace pylith {
  namespace meshio {
    class TestDataWriterHDF5PointsTri3;
    class TestDataWriterHDF5PointsQuad4;
    class TestDataWriterHDF5PointsTet4;
    class TestDataWriterHDF5PointsHex8;
  } // meshio
} // pylith


// ----------------------------------------------------------------------
/// C++ unit testing for DataWriterHDF5
class pylith::meshio::TestDataWriterHDF5PointsTri3 : public TestDataWriterHDF5Points
{ // class TestDataWriterHDF5PointsTri3

  // CPPUNIT TEST SUITE /////////////////////////////////////////////////
  CPPUNIT_TEST_SUITE( TestDataWriterHDF5PointsTri3 );

  CPPUNIT_TEST( testTimeStep );
  CPPUNIT_TEST( testWriteVertexField );

  CPPUNIT_TEST_SUITE_END();

  // PUBLIC METHODS /////////////////////////////////////////////////////
public :

  /// Setup testing data.
  void setUp(void);

}; // class TestDataWriterHDF5PointsTri3


// ----------------------------------------------------------------------
/// C++ unit testing for DataWriterHDF5
class pylith::meshio::TestDataWriterHDF5PointsQuad4 : public TestDataWriterHDF5Points
{ // class TestDataWriterHDF5PointsQuad4

  // CPPUNIT TEST SUITE /////////////////////////////////////////////////
  CPPUNIT_TEST_SUITE( TestDataWriterHDF5PointsQuad4 );

  CPPUNIT_TEST( testTimeStep );
  CPPUNIT_TEST( testWriteVertexField );

  CPPUNIT_TEST_SUITE_END();

  // PUBLIC METHODS /////////////////////////////////////////////////////
public :

  /// Setup testing data.
  void setUp(void);

}; // class TestDataWriterHDF5PointsQuad4


// ----------------------------------------------------------------------
/// C++ unit testing for DataWriterHDF5
class pylith::meshio::TestDataWriterHDF5PointsTet4 : public TestDataWriterHDF5Points
{ // class TestDataWriterHDF5PointsTet4

  // CPPUNIT TEST SUITE /////////////////////////////////////////////////
  CPPUNIT_TEST_SUITE( TestDataWriterHDF5PointsTet4 );

  CPPUNIT_TEST( testTimeStep );
  CPPUNIT_TEST( testWriteVertexField );

  CPPUNIT_TEST_SUITE_END();

  // PUBLIC METHODS /////////////////////////////////////////////////////
public :

  /// Setup testing data.
  void setUp(void);

}; // class TestDataWriterHDF5PointsTet4


// ----------------------------------------------------------------------
/// C++ unit testing for DataWriterHDF5
class pylith::meshio::TestDataWriterHDF5PointsHex8 : public TestDataWriterHDF5Points
{ // class TestDataWriterHDF5PointsHex8

  // CPPUNIT TEST SUITE /////////////////////////////////////////////////
  CPPUNIT_TEST_SUITE( TestDataWriterHDF5PointsHex8 );

  CPPUNIT_TEST( testTimeStep );
  CPPUNIT_TEST( testWriteVertexField );

  CPPUNIT_TEST_SUITE_END();

  // PUBLIC METHODS /////////////////////////////////////////////////////
public :

  /// Setup testing data.
  void setUp(void);

}; // class TestDataWriterHDF5PointsHex8


#endif // pylith_meshio_testdatawriterhdf5pointscases_hh


// End of file 
