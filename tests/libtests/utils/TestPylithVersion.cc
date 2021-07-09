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

#include <portinfo>

#include "TestPylithVersion.hh" // Implementation of class methods

#include "pylith/utils/PylithVersion.hh" // USES PylithVersion

#include <string> // USES std::string()
#include <string.h> // USES strlen()

// ----------------------------------------------------------------------
CPPUNIT_TEST_SUITE_REGISTRATION( pylith::utils::TestPylithVersion );

// ----------------------------------------------------------------------
// Test isRelease()
void
pylith::utils::TestPylithVersion::testIsRelease(void)
{ // testIsRelease
#if PYLITH_RELEASE_VERSION
  CPPUNIT_ASSERT(PylithVersion::isRelease());
#else
  CPPUNIT_ASSERT(!PylithVersion::isRelease());
#endif
} // testIsRelease

// ----------------------------------------------------------------------
// Test version()
void
pylith::utils::TestPylithVersion::testVersion(void)
{ // testVersion
  CPPUNIT_ASSERT_EQUAL(std::string(PYLITH_VERSION), std::string(PylithVersion::version()));
} // testVersion

// ----------------------------------------------------------------------
// Test gitRevision()
void
pylith::utils::TestPylithVersion::testGitRevision(void)
{ // testGitRevision
#if PYLITH_RELEASE_VERSION
  CPPUNIT_ASSERT_EQUAL(std::string("unknown"), std::string(PylithVersion::gitRevision()));
#else
  // Git revision should be of the form vX.X.X-gXXXXX.
  const char* rev = PylithVersion::gitRevision();
  CPPUNIT_ASSERT_EQUAL('v', rev[0]);
#endif
} // testGitRevision

// ----------------------------------------------------------------------
// Test gitHash()
void
pylith::utils::TestPylithVersion::testGitHash(void)
{ // testGitHash
#if PYLITH_RELEASE_VERSION
  CPPUNIT_ASSERT_EQUAL(std::string("unknown"), std::string(PylithVersion::gitHash()));
#else
  // Git hash should contain lower case and numbers.
  const char* hash = PylithVersion::gitHash();
  const int len = strlen(hash);
  for (int i=0; i < len; ++i) {
    const int value = int(hash[i]);
    CPPUNIT_ASSERT((value >= int('0') && value <= int('9')) || (value >= int('a') && value <= int('z')));
  } // for
#endif
} // testGitHash

// ----------------------------------------------------------------------
// Test gitDate()
void
pylith::utils::TestPylithVersion::testGitDate(void)
{ // testGitDate
#if PYLITH_RELEASE_VERSION
  CPPUNIT_ASSERT_EQUAL(std::string("unknown"), std::string(PylithVersion::gitDate()));
#else
  const char* datetime = PylithVersion::gitDate();
  CPPUNIT_ASSERT(strlen(datetime) > 0);
#endif
} // testGitDate

// ----------------------------------------------------------------------
// Test gitBranch()
void
pylith::utils::TestPylithVersion::testGitBranch(void)
{ // testGitBranch
#if PYLITH_RELEASE_VERSION
  CPPUNIT_ASSERT_EQUAL(std::string("unknown"), std::string(PylithVersion::gitBranch()));
#else
  const char* branch = PylithVersion::gitBranch();
  CPPUNIT_ASSERT(strlen(branch) > 0);
#endif
} // testGitBranch


// End of file 
