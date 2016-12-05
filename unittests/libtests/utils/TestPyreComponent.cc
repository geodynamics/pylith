// -*- C++ -*-
//
// ----------------------------------------------------------------------
//
// Brad T. Aagaard, U.S. Geological Survey
// Charles A. Williams, GNS Science
// Matthew G. Knepley, University of Chicago
//
// This code was developed as part of the Computational Infrastructure
// for Geodynamics (http://geodynamics.org).
//
// Copyright (c) 2010-2016 University of California, Davis
//
// See COPYING for license information.
//
// ----------------------------------------------------------------------
//

#include <portinfo>

#include "TestPyreComponent.hh" // Implementation of class methods

#include "pylith/utils/PyreComponent.hh" // USES PyreComponent

#include "pylith/utils/journals.hh" // USES PYLITH_JOURNAL_*

// ----------------------------------------------------------------------
CPPUNIT_TEST_SUITE_REGISTRATION( pylith::utils::TestPyreComponent );

// ----------------------------------------------------------------------
// Test constructor.
void
pylith::utils::TestPyreComponent::testConstructor(void)
{ // testConstructor
    PyreComponent component;

    CPPUNIT_ASSERT_EQUAL(std::string(""), component._name);
} // testConstructor

// ----------------------------------------------------------------------
// Test name().
void
pylith::utils::TestPyreComponent::testName(void)
{ // testName
    PyreComponent component;
    CPPUNIT_ASSERT_EQUAL(std::string(""), std::string(component.name()));

    const std::string& name = "my name";
    component.name(name.c_str());
    CPPUNIT_ASSERT_EQUAL(name, std::string(component.name()));
} // testName

// ----------------------------------------------------------------------
// Test identifier().
void
pylith::utils::TestPyreComponent::testIdentifier(void)
{ // testIdentifier
    PyreComponent component;
    component.name("my component");
    CPPUNIT_ASSERT_EQUAL(std::string("unknown"), std::string(component._identifier));

    const std::string& identifier = "my identifier";
    component.identifier(identifier.c_str());
    CPPUNIT_ASSERT_EQUAL(identifier, std::string(component.identifier()));
} // testIdentifier

// ----------------------------------------------------------------------
// Test PYLITH_JOURNAL_*.
namespace pylith {
    namespace utils {
        class TestJournals : public PyreComponent {
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
pylith::utils::TestPyreComponent::testJournals(void)
{ // testJournals
    const char* name = "test";
    journal::info_t info(name); info.activate();
    journal::debug_t debug(name); debug.activate();
    journal::warning_t warning(name); warning.activate();
    journal::error_t error(name); error.activate();

    TestJournals journals;
    journals.name("test");
    journals.identifier("TestJournals");
    journals.test();
} // testJournals


// End of file
