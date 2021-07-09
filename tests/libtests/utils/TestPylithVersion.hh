// -*- C++ -*-
//
// ----------------------------------------------------------------------
//
// Brad T. Aagaard, U.S. Geological Survey
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

/** @file tests/libtests/utils/TestPylithVersion.hh
 *
 * @brief C++ TestPylithVersion object
 *
 * C++ unit testing for PylithVersion.
 */

#if !defined(pylith_testpylithversion_hh)
#define pylith_testpylithversion_hh

#include <cppunit/extensions/HelperMacros.h>

/// Namespace for pylith package
namespace pylith {
  namespace utils {
    class TestPylithVersion;
  } // utils
} // pylith

/// C++ unit testing for PylithVersion
class pylith::utils::TestPylithVersion : public CppUnit::TestFixture
{ // class TestPylithVersion

  // CPPUNIT TEST SUITE /////////////////////////////////////////////////
  CPPUNIT_TEST_SUITE( TestPylithVersion );

  CPPUNIT_TEST( testIsRelease );
  CPPUNIT_TEST( testVersion );
  CPPUNIT_TEST( testGitRevision );
  CPPUNIT_TEST( testGitHash );
  CPPUNIT_TEST( testGitDate );
  CPPUNIT_TEST( testGitBranch );

  CPPUNIT_TEST_SUITE_END();

  // PUBLIC METHODS /////////////////////////////////////////////////////
public :

  /// Test isRelease()
  void testIsRelease(void);

  /// Test version()
  void testVersion(void);

  /// Test gitRevision()
  void testGitRevision(void);

  /// Test gitHash()
  void testGitHash(void);

  /// Test gitDate()
  void testGitDate(void);

  /// Test gitBranch()
  void testGitBranch(void);

}; // class TestPylithVersion

#endif // pylith_utils_testpylithversion_hh

// End of file 
