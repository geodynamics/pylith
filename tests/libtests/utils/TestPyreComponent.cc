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
// Copyright (c) 2010-2021 University of California, Davis
//
// See LICENSE.md for license information.
//
// ----------------------------------------------------------------------
//

#include <portinfo>

#include "TestPyreComponent.hh" // Implementation of class methods

#include "pylith/utils/PyreComponent.hh" // USES PyreComponent

#include "pylith/utils/journals.hh" // USES PYLITH_COMPONENT_*

// ----------------------------------------------------------------------
CPPUNIT_TEST_SUITE_REGISTRATION(pylith::utils::TestPyreComponent);

// ----------------------------------------------------------------------
// Test constructor.
void
pylith::utils::TestPyreComponent::testConstructor(void) {
    PyreComponent component;

    CPPUNIT_ASSERT_EQUAL(std::string(""), component._name);
} // testConstructor


// ----------------------------------------------------------------------
// Test name().
void
pylith::utils::TestPyreComponent::testName(void) {
    PyreComponent component;
    CPPUNIT_ASSERT_EQUAL(std::string(""), std::string(component.getName()));

    const std::string& name = "my name";
    component.setName(name.c_str());
    CPPUNIT_ASSERT_EQUAL(name, std::string(component.getName()));
} // testName


// ----------------------------------------------------------------------
// Test identifier().
void
pylith::utils::TestPyreComponent::testIdentifier(void) {
    PyreComponent component;
    component.setName("my component");
    CPPUNIT_ASSERT_EQUAL(std::string("unknown"), std::string(component._identifier));

    const std::string& identifier = "my identifier";
    component.setIdentifier(identifier.c_str());
    CPPUNIT_ASSERT_EQUAL(identifier, std::string(component.getIdentifier()));
} // testIdentifier


// ----------------------------------------------------------------------
// Test PYLITH_JOURNAL_*.
namespace pylith {
    namespace utils {
        class TestComponentJournals : public PyreComponent {
public:

            void test(void) {
                PYLITH_COMPONENT_DEBUG("CORRECT: This is a debug message.");
                PYLITH_COMPONENT_INFO("CORRECT: This is an info mesasge.");
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
