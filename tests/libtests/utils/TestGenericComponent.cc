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

#include "TestGenericComponent.hh" // Implementation of class methods

#include "pylith/utils/GenericComponent.hh" // USES GenericComponent

#include "pylith/utils/journals.hh" // USES PYLITH_JOURNAL_*

// ----------------------------------------------------------------------
CPPUNIT_TEST_SUITE_REGISTRATION(pylith::utils::TestGenericComponent);

// ----------------------------------------------------------------------
// Test constructor.
void
pylith::utils::TestGenericComponent::testConstructor(void) {
    GenericComponent component;

    CPPUNIT_ASSERT_EQUAL(std::string(""), component._name);
} // testConstructor


// ----------------------------------------------------------------------
// Test name().
void
pylith::utils::TestGenericComponent::testName(void) {
    GenericComponent component;
    CPPUNIT_ASSERT_EQUAL(std::string(""), std::string(component.getName()));

    const std::string& name = "my name";
    component.setName(name.c_str());
    CPPUNIT_ASSERT_EQUAL(name, std::string(component.getName()));
} // testName


// ----------------------------------------------------------------------
// Test PYLITH_JOURNAL_*.
namespace pylith {
    namespace utils {
        class TestGenericJournals : public GenericComponent {
public:

            void test(void) {
                PYLITH_JOURNAL_DEBUG("CORRECT: This is a debug message.");
                PYLITH_JOURNAL_INFO("CORRECT: This is an info mesasge.");
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
