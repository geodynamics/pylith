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

#include "pylith/utils/PetscVersion.hh" // Test subject

#include "petsc.h" // USES PETSC_VERSION_*

#include <string> // USES std::string()
#include <string.h> // USES strlen()
#include <stdio.h> // USES snprintf()

#include "catch2/catch_test_macros.hpp"

// ------------------------------------------------------------------------------------------------
namespace pylith {
    namespace utils {
        class TestPetscVersion;
    }
}

class pylith::utils::TestPetscVersion {
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

    /// Test gitDate()
    static
    void testGitDate(void);

    /// Test gitBranch()
    static
    void testGitBranch(void);

    /// Test petscDir()
    static
    void testPetscDir(void);

    /// Test petscArch()
    static
    void testPetscArch(void);

}; // class TestPetscVersion

// ------------------------------------------------------------------------------------------------
TEST_CASE("TestPetscVersion::testIsRelease", "[TestPetscVersion]") {
    pylith::utils::TestPetscVersion::testIsRelease();
}
TEST_CASE("TestPetscVersion::testVersion", "[TestPetscVersion]") {
    pylith::utils::TestPetscVersion::testVersion();
}
TEST_CASE("TestPetscVersion::testGitRevision", "[TestPetscVersion]") {
    pylith::utils::TestPetscVersion::testGitRevision();
}
TEST_CASE("TestPetscVersion::testGitDate", "[TestPetscVersion]") {
    pylith::utils::TestPetscVersion::testGitDate();
}
TEST_CASE("TestPetscVersion::testGitBranch", "[TestPetscVersion]") {
    pylith::utils::TestPetscVersion::testGitBranch();
}
TEST_CASE("TestPetscVersion::testPetscDir", "[TestPetscVersion]") {
    pylith::utils::TestPetscVersion::testPetscDir();
}
TEST_CASE("TestPetscVersion::testPetscArch", "[TestPetscVersion]") {
    pylith::utils::TestPetscVersion::testPetscArch();
}

// ------------------------------------------------------------------------------------------------
// Test isRelease()
void
pylith::utils::TestPetscVersion::testIsRelease(void) {
#if PETSC_VERSION_RELEASE
    CHECK(PetscVersion::isRelease());
#else
    CHECK(!PetscVersion::isRelease());
#endif
} // testIsRelease


// ------------------------------------------------------------------------------------------------
// Test version()
void
pylith::utils::TestPetscVersion::testVersion(void) {
    const int maxsize = 64;
    char value[maxsize];
    snprintf(value, maxsize-1, "%d.%d.%d", PETSC_VERSION_MAJOR, PETSC_VERSION_MINOR, PETSC_VERSION_SUBMINOR);
    CHECK(std::string(value) == std::string(PetscVersion::version()));
} // testVersion


// ------------------------------------------------------------------------------------------------
// Test gitRevision()
void
pylith::utils::TestPetscVersion::testGitRevision(void) {
#if PETSC_VERSION_RELEASE
#else
    // Git revision should be of the form vX.X.X-gXXXXX.
    const char* rev = PetscVersion::gitRevision();
    CHECK('v' == rev[0]);
#endif
} // testGitRevision


// ------------------------------------------------------------------------------------------------
// Test gitDate()
void
pylith::utils::TestPetscVersion::testGitDate(void) {
    const char* datetime = PetscVersion::gitDate();
    CHECK(strlen(datetime) > 0);
} // testGitDate


// ------------------------------------------------------------------------------------------------
// Test gitBranch()
void
pylith::utils::TestPetscVersion::testGitBranch(void) {
#if PETSC_VERSION_RELEASE
#else
    const char* branch = PetscVersion::gitBranch();
    CHECK(strlen(branch) > 0);
#endif
} // testGitBranch


// ------------------------------------------------------------------------------------------------
// Test petscDir()
void
pylith::utils::TestPetscVersion::testPetscDir(void) {
    CHECK(strlen(PetscVersion::petscDir()) > 0);
} // testPetscDir


// ------------------------------------------------------------------------------------------------
// Test petscArch()
void
pylith::utils::TestPetscVersion::testPetscArch(void) {
    // Not defined for PETSc prefix installs, so use minimal test.
    CHECK(strlen(PetscVersion::petscArch()) >= 0);
} // testPetscArch


// End of file
