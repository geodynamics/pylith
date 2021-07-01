// -*- C++ -*-
//
// ----------------------------------------------------------------------
//
// Brad T. Aagaard, U.S. Geological Survey
// Charles A. Williams, GNS Science
// Matthew G. Knepley, University at Buffalo
//
// This code was developed as part of the Computational Infrastructure
// for Geodynamics (http://geodynamics.org).
//
// Copyright (c) 2010-2021 University of California, Davis
//
// See LICENSE.md for license information.
//
// ----------------------------------------------------------------------
//

/**
 * @file tests/libtests/meshio/TestHDF5.hh
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
  CPPUNIT_TEST( testHasGroup );
  CPPUNIT_TEST( testHasDataset );
  CPPUNIT_TEST( testGetDatasetDims );
  CPPUNIT_TEST( testGetGroupDatasets );
  CPPUNIT_TEST( testCreateGroup );
  CPPUNIT_TEST( testAttributeScalar );
  CPPUNIT_TEST( testCreateDataset );
  CPPUNIT_TEST( testDatasetChunk );
  CPPUNIT_TEST( testDatasetRawExternal );

  CPPUNIT_TEST( testAttributeString );
  CPPUNIT_TEST( testDatasetString );

  CPPUNIT_TEST_SUITE_END();

  // PUBLIC METHODS /////////////////////////////////////////////////////
public :

  /// Test constructor.
  void testConstructor(void);

  /// Test open() and close().
  void testOpenClose(void);

  /// Test hasGroup().
  void testHasGroup(void);

  /// Test hasDataset().
  void testHasDataset(void);

  /// Test getDatasetDims().
  void testGetDatasetDims(void);

  /// Test getGroupDatasets().
  void testGetGroupDatasets(void);

  /// Test createGroup()
  void testCreateGroup(void);

  /// Test writeAttribute(scalar) and readAttribute(scalar).
  void testAttributeScalar(void);

  /// Test createDataset().
  void testCreateDataset(void);

  /// Test writeDatasetChunk() and readDatasetChunk().
  void testDatasetChunk(void);

  /// Test createDatasetRawExternal() and updateDatasetRawExternal().
  void testDatasetRawExternal(void);

  /// Test writeAttribute(string) and readAttribute(string).
  void testAttributeString(void);

  /// Test writeDataset(string) and readDataset(string).
  void testDatasetString(void);

}; // class TestHDF5

#endif // pylith_meshio_testhdf5_hh

// End of file 
