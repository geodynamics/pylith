// =================================================================================================
// This code is part of PyLith, developed through the Computational Infrastructure
// for Geodynamics (https://github.com/geodynamics/pylith).
//
// Copyright (c) 2010-2023, University of California, Davis and the PyLith Development Team.
// All rights reserved.
//
// See https://mit-license.org/ and LICENSE.md and for license information.
// =================================================================================================

#include <portinfo>

#include "pylith/utils/GenericComponent.hh" // Test subject

#include "pylith/utils/journals.hh" // USES PYLITH_JOURNAL_*

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
// Test PYLITH_JOURNAL_*.
namespace pylith {
    namespace utils {
        class TestGenericJournals : public GenericComponent {
public:

            void test(void) {
                PYLITH_JOURNAL_DEBUG("CORRECT: This is a debug message.");
                PYLITH_JOURNAL_INFO("CORRECT: This is an info mesasge.");
                PYLITH_JOURNAL_INFO_ROOT("CORRECT: This is an info mesasge.");
                PYLITH_JOURNAL_WARNING("CORRECT: This is a warning mesasge.");
                PYLITH_JOURNAL_ERROR("CORRECT: This is an error mesage.");
            } // testJournals

        };
    } // utils
} // pylith

void
pylith::utils::TestGenericComponent::testJournals(void) {
    const char* name = "test";
    pythia::journal::info_t info(name);info.activate();
    pythia::journal::debug_t debug(name);debug.activate();
    pythia::journal::warning_t warning(name);warning.activate();
    pythia::journal::error_t error(name);error.activate();

    TestGenericJournals journals;
    journals.setName("test");
    journals.test();
} // testJournals


// End of file
