// -*- C++ -*-
//
// ----------------------------------------------------------------------
//
// Brad T. Aagaard, U.S. Geological Survey
// Charles A. Williams, GNS Science
// Matthew G. Knepley, University at Buffalo
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

#include "pylith/utils/PyreComponent.hh" // Test subject

#include "pylith/utils/journals.hh" // USES PYLITH_COMPONENT_*

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

    /// Test PYLITH_JOURNAL_DEBUG(), PYLITH_JOURNAL_INFO(), PYLITH_JOURNAL_ERROR().
    static
    void testJournals(void);

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
TEST_CASE("TestPyreComponent::testJournals", "[TestPyreComponent]") {
    pylith::utils::TestPyreComponent::testJournals();
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
                PYLITH_COMPONENT_DEBUG("CORRECT: This is a debug message.");
                PYLITH_COMPONENT_INFO("CORRECT: This is an info mesasge.");
                PYLITH_COMPONENT_INFO_ROOT("CORRECT: This is an info mesasge.");
                PYLITH_COMPONENT_WARNING("CORRECT: This is a warning mesasge.");
                PYLITH_COMPONENT_ERROR("CORRECT: This is an error mesage.");
            } // testJournals

        };
    } // utils
} // pylith

void
pylith::utils::TestPyreComponent::testJournals(void) {
    const char* name = "test";
    pythia::journal::info_t info(name);info.activate();
    pythia::journal::debug_t debug(name);debug.activate();
    pythia::journal::warning_t warning(name);warning.activate();
    pythia::journal::error_t error(name);error.activate();

    TestComponentJournals journals;
    journals.setName("test");
    journals.setIdentifier("TestJournals");
    journals.test();
} // testJournals


// End of file
