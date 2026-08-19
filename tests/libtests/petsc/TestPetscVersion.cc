// =================================================================================================
// This code is part of PyLith, developed through the Computational Infrastructure
// for Geodynamics (https://github.com/geodynamics/pylith).
//
// Copyright (c) 2010-2026, University of California, Davis and the PyLith Development Team.
// All rights reserved.
//
// See https://mit-license.org/ and LICENSE.md and for license information.
// =================================================================================================
//

#include <portinfo>

#include "pylith/petsc/PetscVersion.hh" // Test subject

#include "petsc.h" // USES PETSC_VERSION_*

#include <string> // USES std::string()
#include <string.h> // USES strlen()
#include <stdio.h> // USES snprintf()

#include "catch2/catch_test_macros.hpp"

// ------------------------------------------------------------------------------------------------
namespace pylith {
    namespace petsc {
        class TestPetscVersion;
    }
}

class pylith::petsc::TestPetscVersion {
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
    pylith::petsc::TestPetscVersion::testIsRelease();
}
TEST_CASE("TestPetscVersion::testVersion", "[TestPetscVersion]") {
    pylith::petsc::TestPetscVersion::testVersion();
}
TEST_CASE("TestPetscVersion::testGitRevision", "[TestPetscVersion]") {
    pylith::petsc::TestPetscVersion::testGitRevision();
}
TEST_CASE("TestPetscVersion::testGitDate", "[TestPetscVersion]") {
    pylith::petsc::TestPetscVersion::testGitDate();
}
TEST_CASE("TestPetscVersion::testGitBranch", "[TestPetscVersion]") {
    pylith::petsc::TestPetscVersion::testGitBranch();
}
TEST_CASE("TestPetscVersion::testPetscDir", "[TestPetscVersion]") {
    pylith::petsc::TestPetscVersion::testPetscDir();
}
TEST_CASE("TestPetscVersion::testPetscArch", "[TestPetscVersion]") {
    pylith::petsc::TestPetscVersion::testPetscArch();
}

// ------------------------------------------------------------------------------------------------
// Test isRelease()
void
pylith::petsc::TestPetscVersion::testIsRelease(void) {
#if PETSC_VERSION_RELEASE
    CHECK(PetscVersion::isRelease());
#else
    CHECK(!PetscVersion::isRelease());
#endif
} // testIsRelease


// ------------------------------------------------------------------------------------------------
// Test version()
void
pylith::petsc::TestPetscVersion::testVersion(void) {
    const int maxsize = 64;
    char value[maxsize];
    snprintf(value, maxsize-1, "%d.%d.%d", PETSC_VERSION_MAJOR, PETSC_VERSION_MINOR, PETSC_VERSION_SUBMINOR);
    CHECK(std::string(value) == std::string(PetscVersion::version()));
} // testVersion


// ------------------------------------------------------------------------------------------------
// Test gitRevision()
void
pylith::petsc::TestPetscVersion::testGitRevision(void) {
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
pylith::petsc::TestPetscVersion::testGitDate(void) {
    const char* datetime = PetscVersion::gitDate();
    CHECK(strlen(datetime) > 0);
} // testGitDate


// ------------------------------------------------------------------------------------------------
// Test gitBranch()
void
pylith::petsc::TestPetscVersion::testGitBranch(void) {
#if PETSC_VERSION_RELEASE
#else
    const char* branch = PetscVersion::gitBranch();
    CHECK(strlen(branch) > 0);
#endif
} // testGitBranch


// ------------------------------------------------------------------------------------------------
// Test petscDir()
void
pylith::petsc::TestPetscVersion::testPetscDir(void) {
    CHECK(strlen(PetscVersion::petscDir()) > 0);
} // testPetscDir


// ------------------------------------------------------------------------------------------------
// Test petscArch()
void
pylith::petsc::TestPetscVersion::testPetscArch(void) {
    // Not defined for PETSc prefix installs, so empty string is valid.
} // testPetscArch


// End of file
