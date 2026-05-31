// =================================================================================================
// This code is part of PyLith, developed through the Computational Infrastructure
// for Geodynamics (https://github.com/geodynamics/pylith).
//
// Copyright (c) 2010-2025, University of California, Davis and the PyLith Development Team.
// All rights reserved.
//
// See https://mit-license.org/ and LICENSE.md and for license information.
// =================================================================================================

#include <portinfo>

#include "pylith/utils/GenericComponent.hh" // Test subject

#include "pylith/utils/journals.hh" // USES PYLITH_JOURNAL_*
#include "pylith/utils/Exceptions.hh" // USES Exception

#include "catch2/catch_test_macros.hpp"

// ------------------------------------------------------------------------------------------------
namespace pylith {
    namespace utils {
        class TestGenericComponent;
    }
}

class pylith::utils::TestGenericComponent {
    // PUBLIC METHODS /////////////////////////////////////////////////////////////////////////////
public:

    /// Test constructor.
    static
    void testConstructor(void);

    /// Test name().
    static
    void testName(void);

    /// Test PYLITH_JOURNAL_DEBUG(), PYLITH_JOURNAL_INFO(), PYLITH_JOURNAL_ERROR().
    static
    void testJournals(void);

}; // class TestGenericComponent

// ------------------------------------------------------------------------------------------------
TEST_CASE("TestGenericComponent::testConstructor", "[TestGenericComponent]") {
    pylith::utils::TestGenericComponent::testConstructor();
}
TEST_CASE("TestGenericComponent::testName", "[TestGenericComponent]") {
    pylith::utils::TestGenericComponent::testName();
}
TEST_CASE("TestGenericComponent::testJournals", "[TestGenericComponent]") {
    pylith::utils::TestGenericComponent::testJournals();
}

// ------------------------------------------------------------------------------------------------
// Test constructor.
void
pylith::utils::TestGenericComponent::testConstructor(void) {
    GenericComponent component;

    CHECK(std::string("") == component._name);
} // testConstructor


// ------------------------------------------------------------------------------------------------
// Test name().
void
pylith::utils::TestGenericComponent::testName(void) {
    GenericComponent component;
    CHECK(std::string("") == std::string(component.getName()));

    const std::string& name = "my name";
    component.setName(name.c_str());
    CHECK(name == std::string(component.getName()));
} // testName


// ------------------------------------------------------------------------------------------------
void
pylith::utils::TestGenericComponent::testJournals(void) {
    pythia::journal::info_t info(pylith::journal::application_flow);info.activate();
    pythia::journal::info_t info2(pylith::journal::about);info2.activate();
    pythia::journal::debug_t debug(pylith::journal::mms_test);debug.activate();
    pythia::journal::warning_t warning(pylith::journal::deprecation);warning.activate();
    pythia::journal::error_t error(pylith::journal::user_input);error.activate();
    pythia::journal::error_t error2(pylith::journal::logic);error2.activate();

    PYLITH_INFO_ROOT(pylith::journal::application_flow, "CORRECT: This is an info message.");
    PYLITH_INFO(pylith::journal::about, "CORRECT: This is an info message.");
    PYLITH_DEBUG(pylith::journal::mms_test, "CORRECT: This is a debug message.");
    PYLITH_WARNING(pylith::journal::deprecation, "CORRECT: This is a warning message.");

    // CHECK_THROWS_AS() does not seem to work with PYLITH_COMPONENT* macros.
    try {
        PYLITH_ERROR(pylith::IOError, pylith::journal::user_input, "CORRECT: This is an error message.");
    } catch (const pylith::IOError& err) {}
    try {
        PYLITH_FIREWALL(pylith::InternalLogicError, pylith::journal::logic, "CORRECT: This is an error message.");
    } catch (const pylith::InternalLogicError& err) {}
} // testJournals


// End of file
