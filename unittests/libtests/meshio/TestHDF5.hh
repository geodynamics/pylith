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
 * @file unittests/libtests/meshio/TestHDF5.hh
 *
 * @brief C++ TestHDF5 object
 *
 * C++ unit testing for HDF5.
 */

#if !defined(pylith_meshio_testhdf5_hh)
#define pylith_meshio_testhdf5_hh

#include <cppunit/extensions/HelperMacros.h>

/// Namespace for pylith package
namespace pylith {
  namespace meshio {
    class TestHDF5;
  } // meshio
} // pylith

/// C++ unit testing for HDF5
class pylith::meshio::TestHDF5 : public CppUnit::TestFixture
{ // class TestHDF5

  // CPPUNIT TEST SUITE /////////////////////////////////////////////////
  CPPUNIT_TEST_SUITE( TestHDF5 );

  CPPUNIT_TEST( testConstructor );
  CPPUNIT_TEST( testOpenClose );
  CPPUNIT_TEST( testCreateGroup );
  CPPUNIT_TEST( testAttributeScalar );
  CPPUNIT_TEST( testAttributeString );
  CPPUNIT_TEST( testCreateDataset );
  CPPUNIT_TEST( testCreateDatasetRawExternal );
  CPPUNIT_TEST( testDatasetChunk );

  CPPUNIT_TEST_SUITE_END();

  // PUBLIC METHODS /////////////////////////////////////////////////////
public :

  /// Test constructor.
  void testConstructor(void);

  /// Test open() and close().
  void testOpenClose(void);

  /// Test createGroup()
  void testCreateGroup(void);

  /// Test writeAttribute(scalar) and readAttribute(scalar).
  void testAttributeScalar(void);

  /// Test writeAttribute(string) and readAttribute(string).
  void testAttributeString(void);

  /// Test createDataset().
  void testCreateDataset(void);

  /// Test createDatasetRawExternal().
  void testCreateDatasetRawExternal(void);

  /// Test writeDatasetChunk() and readDatasetChunk().
  void testDatasetChunk(void);

}; // class TestHDF5

#endif // pylith_meshio_testhdf5_hh

// End of file 
