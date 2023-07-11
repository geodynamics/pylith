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

#include "pylith/utils/DependenciesVersion.hh" // Test subject

#include <string.h> // USES strlen()

#include "catch2/catch_test_macros.hpp"

// ------------------------------------------------------------------------------------------------
namespace pylith {
    namespace utils {
        class TestDependenciesVersion;
    }
}

class pylith::utils::TestDependenciesVersion {
    // PUBLIC METHODS /////////////////////////////////////////////////////////////////////////////
public:

    /// Test mpiVersion()
    static
    void testMPIVersion(void);

    /// Test mpiImplementation()
    static
    void testMPIImplementation(void);

    /// Test mpiStandard()
    static
    void testMPIStandard(void);

    /// Test netcdfVersion()
    static
    void testNetCDFVersion(void);

    /// Test hdf5Version()
    static
    void testHDF5Version(void);

}; // class TestDependenciesVersion

// ------------------------------------------------------------------------------------------------
TEST_CASE("TestDependenciesVersion::testMPIVersion", "[TestDependenciesVersion]") {
    pylith::utils::TestDependenciesVersion::testMPIVersion();
}
TEST_CASE("TestDependenciesVersion::testMPIImplementation", "[TestDependenciesVersion]") {
    pylith::utils::TestDependenciesVersion::testMPIImplementation();
}
TEST_CASE("TestDependenciesVersion::testMPIStandard", "[TestDependenciesVersion]") {
    pylith::utils::TestDependenciesVersion::testMPIStandard();
}
TEST_CASE("TestDependenciesVersion::testNetCDFVersion", "[TestDependenciesVersion]") {
    pylith::utils::TestDependenciesVersion::testNetCDFVersion();
}
TEST_CASE("TestDependenciesVersion::testHDF5Version", "[TestDependenciesVersion]") {
    pylith::utils::TestDependenciesVersion::testHDF5Version();
}

// ------------------------------------------------------------------------------------------------
// Test mpiVersion()
void
pylith::utils::TestDependenciesVersion::testMPIVersion(void) {
    const char* v = DependenciesVersion::mpiVersion();
    CHECK(strlen(v) > 0);
} // testMPIVersion


// ------------------------------------------------------------------------------------------------
// Test mpiVersion()
void
pylith::utils::TestDependenciesVersion::testMPIImplementation(void) {
    const char* v = DependenciesVersion::mpiImplementation();
    CHECK(strlen(v) > 0);
} // testMPIImplementation


// ------------------------------------------------------------------------------------------------
// Test mpiStandard()
void
pylith::utils::TestDependenciesVersion::testMPIStandard(void) {
    const char* v = DependenciesVersion::mpiStandard();
    CHECK(strlen(v) > 0);
} // testMPIStandard


// ------------------------------------------------------------------------------------------------
// Test netcdfVersion()
void
pylith::utils::TestDependenciesVersion::testNetCDFVersion(void) {
    const char* v = DependenciesVersion::netcdfVersion();
    CHECK(strlen(v) > 0);
} // testNetCDFVersion


// ------------------------------------------------------------------------------------------------
// Test hdf5Version()
void
pylith::utils::TestDependenciesVersion::testHDF5Version(void) {
    const char* v = DependenciesVersion::hdf5Version();
    CHECK(strlen(v) > 0);
} // testHDF5Version


// End of file
