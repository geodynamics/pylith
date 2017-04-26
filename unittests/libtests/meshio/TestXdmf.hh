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
// Copyright (c) 2010-2017 University of California, Davis
//
// See COPYING for license information.
//
// ----------------------------------------------------------------------
//

/**
 * @file unittests/libtests/meshio/TestXdmf.hh
 *
 * @brief C++ TestXdmf object
 *
 * C++ unit testing for Xdmf.
 */

#if !defined(pylith_meshio_testxdmf_hh)
#define pylith_meshio_testxdmf_hh

#include <cppunit/extensions/HelperMacros.h>

/// Namespace for pylith package
namespace pylith {
  namespace meshio {
    class TestXdmf;

    class TestXdmf_Data;
  } // meshio
} // pylith

// ======================================================================
/// C++ unit testing for Xdmf
class pylith::meshio::TestXdmf : public CppUnit::TestFixture {

  // CPPUNIT TEST SUITE /////////////////////////////////////////////////
  CPPUNIT_TEST_SUITE( TestXdmf );

  CPPUNIT_TEST( testConstructor );
  CPPUNIT_TEST( testWrite );

  CPPUNIT_TEST_SUITE_END();

  // PUBLIC METHODS /////////////////////////////////////////////////////
public :

  /// Setup testing data.
  void setUp(void);

  /// Deallocate testing data.
  void tearDown(void);

  /// Test constructor.
  void testConstructor(void);

  /// Test write().
  void testWrite(void);

  // PROTECTED MEMBDERS /////////////////////////////////////////////////
protected:

  TestXdmf_Data* _data; ///< Data for testing.

}; // class TestXdmf


// ======================================================================
class pylith::meshio::TestXdmf_Data {

    // PUBLIC METHODS ///////////////////////////////////////////////////////
public:

    /// Constructor
    TestXdmf_Data(void);

    /// Destructor
    ~TestXdmf_Data(void);

    // PUBLIC MEMBERS ///////////////////////////////////////////////////////
public:

    const char* filenameHDF5;
    const char* filenameXdmf;

}; // class TestXdmf_Data

#endif // pylith_meshio_testxdmf_hh

// End of file 
