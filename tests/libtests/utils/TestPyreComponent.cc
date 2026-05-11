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

#include "pylith/utils/PyreComponent.hh" // Test subject

#include "pylith/utils/journals.hh" // USES PYLITH_COMPONENT_*
#include "pylith/utils/Exceptions.hh" // USES Exception

#include "catch2/catch_test_macros.hpp"

// ------------------------------------------------------------------------------------------------
namespace pylith {
    namespace utils {
        class TestPyreComponent;
    }
}

class pylith::utils::TestPyreComponent {
    // PUBLIC METHODS /////////////////////////////////////////////////////////////////////////////
public:

    /// Test constructor.
    static
    void testConstructor(void);

    /// Test name().
    static
    void testName(void);

    /// Test identifier().
    static
    void testIdentifier(void);

    /// Test component journal macros
    static
    void testComponentJournalMacros(void);

    /// Test journal macros
    static
    void testJournalMacros(void);

}; // TestPyreComponent

// ------------------------------------------------------------------------------------------------
TEST_CASE("TestPyreComponent::testConstructor", "[TestPyreComponent]") {
    pylith::utils::TestPyreComponent::testConstructor();
}
TEST_CASE("TestPyreComponent::testName", "[TestPyreComponent]") {
    pylith::utils::TestPyreComponent::testName();
}
TEST_CASE("TestPyreComponent::testIdentifier", "[TestPyreComponent]") {
    pylith::utils::TestPyreComponent::testIdentifier();
}
TEST_CASE("TestPyreComponent::testComponentJournalMacros", "[TestPyreComponent]") {
    pylith::utils::TestPyreComponent::testComponentJournalMacros();
}
TEST_CASE("TestPyreComponent::testJournalMacros", "[TestPyreComponent]") {
    pylith::utils::TestPyreComponent::testJournalMacros();
}

// ------------------------------------------------------------------------------------------------
// Test constructor.
void
pylith::utils::TestPyreComponent::testConstructor(void) {
    PyreComponent component;

    CHECK(std::string("") == component._name);
} // testConstructor


// ------------------------------------------------------------------------------------------------
// Test name().
void
pylith::utils::TestPyreComponent::testName(void) {
    PyreComponent component;
    CHECK(std::string("") == std::string(component.getName()));

    const std::string& name = "my name";
    component.setName(name.c_str());
    CHECK(name == std::string(component.getName()));
} // testName


// ------------------------------------------------------------------------------------------------
// Test identifier().
void
pylith::utils::TestPyreComponent::testIdentifier(void) {
    PyreComponent component;
    component.setName("my component");
    CHECK(std::string("unknown") == std::string(component._identifier));

    const std::string& identifier = "my identifier";
    component.setIdentifier(identifier.c_str());
    CHECK(identifier == std::string(component.getIdentifier()));
} // testIdentifier


// ------------------------------------------------------------------------------------------------
// Test PYLITH_JOURNAL_*.
namespace pylith {
    namespace utils {
        class TestComponentJournals : public PyreComponent {
public:

            void test(void) {
                PYLITH_COMPONENT_INFO_ROOT_NEW(pylith::journal::application_flow, "CORRECT: This is an info message.");
                PYLITH_COMPONENT_INFO_NEW(pylith::journal::about, "CORRECT: This is an info message.");
                PYLITH_COMPONENT_DEBUG_NEW(pylith::journal::mms_test, "CORRECT: This is a debug message.");
                PYLITH_COMPONENT_WARNING_NEW(pylith::journal::deprecation, "CORRECT: This is a warning message.");

                // CHECK_THROWS_AS() does not seem to work with PYLITH_COMPONENT* macros.
                try {
                    PYLITH_COMPONENT_ERROR_NEW(pylith::IOError, pylith::journal::input_error, "CORRECT: This is an error message.");
                } catch (const pylith::IOError& err) {}
                try {
                    PYLITH_COMPONENT_FIREWALL(pylith::InternalLogicError, pylith::journal::logic_error, "CORRECT: This is an error message.");
                } catch (const pylith::InternalLogicError& err) {}

            } // testJournals

        };

    } // utils
} // pylith


// ------------------------------------------------------------------------------------------------
void
pylith::utils::TestPyreComponent::testComponentJournalMacros(void) {
    pythia::journal::info_t info(pylith::journal::application_flow);info.activate();
    pythia::journal::info_t info2(pylith::journal::about);info2.activate();
    pythia::journal::debug_t debug(pylith::journal::mms_test);debug.activate();
    pythia::journal::warning_t warning(pylith::journal::deprecation);warning.activate();
    pythia::journal::error_t error(pylith::journal::input_error);error.activate();
    pythia::journal::error_t error2(pylith::journal::logic_error);error2.activate();

    TestComponentJournals journals;
    journals.setName("test");
    journals.setIdentifier("TestJournals");
    journals.test();
} // testComponentJournalMacros


// ------------------------------------------------------------------------------------------------
void
pylith::utils::TestPyreComponent::testJournalMacros(void) {
    PYLITH_INFO_ROOT_NEW(pylith::journal::application_flow, "CORRECT: This is an info message.");
    PYLITH_INFO_NEW(pylith::journal::about, "CORRECT: This is an info message.");
    PYLITH_DEBUG_NEW(pylith::journal::mms_test, "CORRECT: This is a debug message.");
    PYLITH_WARNING_NEW(pylith::journal::deprecation, "CORRECT: This is a warning message.");

    // CHECK_THROWS_AS() does not seem to work with PYLITH_COMPONENT* macros.
    try {
        PYLITH_ERROR_NEW(pylith::IOError, pylith::journal::input_error, "CORRECT: This is an error message.");
    } catch (const pylith::IOError& err) {}
    try {
        PYLITH_FIREWALL(pylith::InternalLogicError, pylith::journal::logic_error, "CORRECT: This is an error message.");
    } catch (const pylith::InternalLogicError& err) {}
} // testJournalMacros


// End of file
