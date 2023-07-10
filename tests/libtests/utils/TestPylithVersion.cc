// -*- C++ -*-
//
// ----------------------------------------------------------------------
//
// Brad T. Aagaard, U.S. Geological Survey
//
// This code was developed as part of the Computational Infrastructure
// for Geodynamics (http://geodynamics.org).
//
// Copyright (c) 2010-2022 University of California, Davis
//
// See LICENSE.md for license information.
//
// ----------------------------------------------------------------------
//

#include <portinfo>

#include "pylith/utils/PylithVersion.hh" // Test subject

#include <string> // USES std::string()
#include <string.h> // USES strlen()

#include "catch2/catch_test_macros.hpp"

// ------------------------------------------------------------------------------------------------
namespace pylith {
    namespace utils {
        class TestPylithVersion;
    }
}

class pylith::utils::TestPylithVersion {
    // PUBLIC METHODS /////////////////////////////////////////////////////////////////////////////
public:

    /// Test isRelease()
    static
    void testIsRelease(void);

    /// Test version()
    static
    void testVersion(void);

    /// Test gitRevision()
    static
    void testGitRevision(void);

    /// Test gitHash()
    static
    void testGitHash(void);

    /// Test gitDate()
    static
    void testGitDate(void);

    /// Test gitBranch()
    static
    void testGitBranch(void);

}; // class TestPylithVersion

// ------------------------------------------------------------------------------------------------
TEST_CASE("TestPylithVersion::testIsRelease", "[TestPylithVersion]") {
    pylith::utils::TestPylithVersion::testIsRelease();
}
TEST_CASE("TestPylithVersion::testVersion", "[TestPylithVersion]") {
    pylith::utils::TestPylithVersion::testVersion();
}
TEST_CASE("TestPylithVersion::testGitRevision", "[TestPylithVersion]") {
    pylith::utils::TestPylithVersion::testGitRevision();
}
TEST_CASE("TestPylithVersion::testGitHash", "[TestPylithVersion]") {
    pylith::utils::TestPylithVersion::testGitHash();
}
TEST_CASE("TestPylithVersion::testGitDate", "[TestPylithVersion]") {
    pylith::utils::TestPylithVersion::testGitDate();
}
TEST_CASE("TestPylithVersion::testGitBranch", "[TestPylithVersion]") {
    pylith::utils::TestPylithVersion::testGitBranch();
}

// ------------------------------------------------------------------------------------------------
// Test isRelease()
void
pylith::utils::TestPylithVersion::testIsRelease(void) { // testIsRelease
#if PYLITH_RELEASE_VERSION
    CHECK(PylithVersion::isRelease());
#else
    CHECK(!PylithVersion::isRelease());
#endif
} // testIsRelease


// ------------------------------------------------------------------------------------------------
// Test version()
void
pylith::utils::TestPylithVersion::testVersion(void) { // testVersion
    CHECK(std::string(PYLITH_VERSION) == std::string(PylithVersion::version()));
} // testVersion


// ------------------------------------------------------------------------------------------------
// Test gitRevision()
void
pylith::utils::TestPylithVersion::testGitRevision(void) { // testGitRevision
#if PYLITH_RELEASE_VERSION
    CHECK(std::string("unknown") == std::string(PylithVersion::gitRevision()));
#else
    // Git revision should be of the form vX.X.X-gXXXXX.
    const char* rev = PylithVersion::gitRevision();
    CHECK('v' == rev[0]);
#endif
} // testGitRevision


// ------------------------------------------------------------------------------------------------
// Test gitHash()
void
pylith::utils::TestPylithVersion::testGitHash(void) { // testGitHash
#if PYLITH_RELEASE_VERSION
    CHECK(std::string("unknown") == std::string(PylithVersion::gitHash()));
#else
    // Git hash should contain lower case and numbers.
    const char* hash = PylithVersion::gitHash();
    const int len = strlen(hash);
    for (int i = 0; i < len; ++i) {
        const int value = int(hash[i]);
        assert((value >= int('0') && value <= int('9')) || (value >= int('a') && value <= int('z')));
    } // for
#endif
} // testGitHash


// ------------------------------------------------------------------------------------------------
// Test gitDate()
void
pylith::utils::TestPylithVersion::testGitDate(void) { // testGitDate
#if PYLITH_RELEASE_VERSION
    CHECK(std::string("unknown") == std::string(PylithVersion::gitDate()));
#else
    const char* datetime = PylithVersion::gitDate();
    CHECK(strlen(datetime) > 0);
#endif
} // testGitDate


// ------------------------------------------------------------------------------------------------
// Test gitBranch()
void
pylith::utils::TestPylithVersion::testGitBranch(void) { // testGitBranch
#if PYLITH_RELEASE_VERSION
    CHECK(std::string("unknown") == std::string(PylithVersion::gitBranch()));
#else
    const char* branch = PylithVersion::gitBranch();
    CHECK(strlen(branch) > 0);
#endif
} // testGitBranch


// End of file
