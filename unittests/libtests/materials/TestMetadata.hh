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
 * @file unittests/libtests/materials/TestMetadata.hh
 *
 * @brief C++ TestMetadata object
 *
 * C++ unit testing for Material.
 */

#if !defined(pylith_materials_testmetadata_hh)
#define pylith_materials_testmetadata_hh

#include <cppunit/extensions/HelperMacros.h>

#include "pylith/materials/materialsfwd.hh" // forward declarations

/// Namespace for pylith package
namespace pylith {
  namespace materials {
    class TestMetadata;
  } // materials
} // pylith

/// C++ unit testing for Material
class pylith::materials::TestMetadata : public CppUnit::TestFixture
{ // class TestMetadata

  // CPPUNIT TEST SUITE /////////////////////////////////////////////////
  CPPUNIT_TEST_SUITE( TestMetadata );

  CPPUNIT_TEST( testConstructor );
  CPPUNIT_TEST( testCopyConstructor );
  CPPUNIT_TEST( testProperties );
  CPPUNIT_TEST( testStateVars );
  CPPUNIT_TEST( testDBProperties );
  CPPUNIT_TEST( testDBStateVars );

  CPPUNIT_TEST_SUITE_END();

  // PUBLIC METHODS /////////////////////////////////////////////////////
public :

  /// Setup test data.
  void setUp(void);

  /// Tear down test data.
  void tearDown(void);

  /// Test constructor.
  void testConstructor(void);

  /// Test copy constructor.
  void testCopyConstructor(void);

  /// Test numProperties() and getProperty().
  void testProperties(void);

  /// Test numStateVars() and getStateVar().
  void testStateVars(void);

  /// Test dbProperties() and numDBProperties().
  void testDBProperties(void);

  /// Test dbStateVars() and numDBStateVars().
  void testDBStateVars(void);

  // PRIVATE MEMBERS ////////////////////////////////////////////////////
private :

  Metadata* _metadata; ///< Object for testing

}; // class TestMetadata

#endif // pylith_materials_testmetadata_hh

// End of file 
