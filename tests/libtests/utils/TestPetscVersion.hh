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

/** @file tests/libtests/utils/TestPetscVersion.hh
 *
 * @brief C++ TestPetscVersion object
 *
 * C++ unit testing for PetscVersion.
 */

#if !defined(pylith_testpetscversion_hh)
#define pylith_testpetscversion_hh

#include <cppunit/extensions/HelperMacros.h>

/// Namespace for pylith package
namespace pylith {
  namespace utils {
    class TestPetscVersion;
  } // utils
} // pylith

/// C++ unit testing for PetscVersion
class pylith::utils::TestPetscVersion : public CppUnit::TestFixture
{ // class TestPetscVersion

  // CPPUNIT TEST SUITE /////////////////////////////////////////////////
  CPPUNIT_TEST_SUITE( TestPetscVersion );

  CPPUNIT_TEST( testIsRelease );
  CPPUNIT_TEST( testVersion );
  CPPUNIT_TEST( testGitRevision );
  CPPUNIT_TEST( testGitDate );
  CPPUNIT_TEST( testGitBranch );
  CPPUNIT_TEST( testPetscDir );
  CPPUNIT_TEST( testPetscArch );

  CPPUNIT_TEST_SUITE_END();

  // PUBLIC METHODS /////////////////////////////////////////////////////
public :

  /// Test isRelease()
  void testIsRelease(void);

  /// Test version()
  void testVersion(void);

  /// Test gitRevision()
  void testGitRevision(void);

  /// Test gitDate()
  void testGitDate(void);

  /// Test gitBranch()
  void testGitBranch(void);

  /// Test petscDir()
  void testPetscDir(void);

  /// Test petscArch()
  void testPetscArch(void);

}; // class TestPetscVersion

#endif // pylith_utils_testpetscversion_hh

// End of file 
