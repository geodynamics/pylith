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
 * @file unittests/libtests/meshio/TestDataWriterVTKPoints.hh
 *
 * @brief C++ TestDataWriterVTKPoints object
 *
 * C++ unit testing for DataWriterVTKPoints.
 */

#if !defined(pylith_meshio_testdatawritervtkpoints_hh)
#define pylith_meshio_testdatawritervtkpoints_hh

#include "TestDataWriterVTK.hh" // ISA TestDataWriterVTK
#include "TestDataWriterPoints.hh" // ISA TestDataWriterPoints

#include "pylith/topology/topologyfwd.hh" // USES Mesh, Field

#include <cppunit/extensions/HelperMacros.h>

/// Namespace for pylith package
namespace pylith {
  namespace meshio {
    class TestDataWriterVTKPoints;
  } // meshio
} // pylith

/// C++ unit testing for DataWriterVTK
class pylith::meshio::TestDataWriterVTKPoints : public TestDataWriterVTK,
						public TestDataWriterPoints,
						public CppUnit::TestFixture
{ // class TestDataWriterVTKPoints

  // CPPUNIT TEST SUITE /////////////////////////////////////////////////
  CPPUNIT_TEST_SUITE( TestDataWriterVTKPoints );

  CPPUNIT_TEST( testConstructor );

  CPPUNIT_TEST_SUITE_END();

  // PUBLIC METHODS /////////////////////////////////////////////////////
public :

  /// Setup testing data.
  void setUp(void);

  /// Tear down testing data.
  void tearDown(void);

  /// Test constructor
  void testConstructor(void);

  /// Test openTimeStep() and closeTimeStep()
  void testTimeStep(void);

  /// Test writeVertexField.
  void testWriteVertexField(void);

}; // class TestDataWriterVTKPoints

#endif // pylith_meshio_testdatawritervtkpoints_hh


// End of file 
