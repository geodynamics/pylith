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
 * @file unittests/libtests/meshio/TestDataWriterVTKPointsCases.hh
 *
 * @brief C++ unit testing for DataWriterVTK output interpolated to
 * points with various mesh cell types.
 */

#if !defined(pylith_meshio_testdatawritervtkpointscases_hh)
#define pylith_meshio_testdatawritervtkpointscases_hh

#include "TestDataWriterVTKPoints.hh"

/// Namespace for pylith package
namespace pylith {
  namespace meshio {
    class TestDataWriterVTKPointsTri3;
    class TestDataWriterVTKPointsQuad4;
    class TestDataWriterVTKPointsTet4;
    class TestDataWriterVTKPointsHex8;
  } // meshio
} // pylith


// ----------------------------------------------------------------------
/// C++ unit testing for DataWriterVTK
class pylith::meshio::TestDataWriterVTKPointsTri3 : public TestDataWriterVTKPoints
{ // class TestDataWriterVTKPointsTri3

  // CPPUNIT TEST SUITE /////////////////////////////////////////////////
  CPPUNIT_TEST_SUITE( TestDataWriterVTKPointsTri3 );

  CPPUNIT_TEST( testTimeStep );
  CPPUNIT_TEST( testWriteVertexField );

  CPPUNIT_TEST_SUITE_END();

  // PUBLIC METHODS /////////////////////////////////////////////////////
public :

  /// Setup testing data.
  void setUp(void);

}; // class TestDataWriterVTKPointsTri3


// ----------------------------------------------------------------------
/// C++ unit testing for DataWriterVTK
class pylith::meshio::TestDataWriterVTKPointsQuad4 : public TestDataWriterVTKPoints
{ // class TestDataWriterVTKPointsQuad4

  // CPPUNIT TEST SUITE /////////////////////////////////////////////////
  CPPUNIT_TEST_SUITE( TestDataWriterVTKPointsQuad4 );

  CPPUNIT_TEST( testTimeStep );
  CPPUNIT_TEST( testWriteVertexField );

  CPPUNIT_TEST_SUITE_END();

  // PUBLIC METHODS /////////////////////////////////////////////////////
public :

  /// Setup testing data.
  void setUp(void);

}; // class TestDataWriterVTKPointsQuad4


// ----------------------------------------------------------------------
/// C++ unit testing for DataWriterVTK
class pylith::meshio::TestDataWriterVTKPointsTet4 : public TestDataWriterVTKPoints
{ // class TestDataWriterVTKPointsTet4

  // CPPUNIT TEST SUITE /////////////////////////////////////////////////
  CPPUNIT_TEST_SUITE( TestDataWriterVTKPointsTet4 );

  CPPUNIT_TEST( testTimeStep );
  CPPUNIT_TEST( testWriteVertexField );

  CPPUNIT_TEST_SUITE_END();

  // PUBLIC METHODS /////////////////////////////////////////////////////
public :

  /// Setup testing data.
  void setUp(void);

}; // class TestDataWriterVTKPointsTet4


// ----------------------------------------------------------------------
/// C++ unit testing for DataWriterVTK
class pylith::meshio::TestDataWriterVTKPointsHex8 : public TestDataWriterVTKPoints
{ // class TestDataWriterVTKPointsHex8

  // CPPUNIT TEST SUITE /////////////////////////////////////////////////
  CPPUNIT_TEST_SUITE( TestDataWriterVTKPointsHex8 );

  CPPUNIT_TEST( testTimeStep );
  CPPUNIT_TEST( testWriteVertexField );

  CPPUNIT_TEST_SUITE_END();

  // PUBLIC METHODS /////////////////////////////////////////////////////
public :

  /// Setup testing data.
  void setUp(void);

}; // class TestDataWriterVTKPointsHex8


#endif // pylith_meshio_testdatawritervtkpointscases_hh


// End of file 
