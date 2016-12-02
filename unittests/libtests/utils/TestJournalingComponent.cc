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

#include "TestJournalingComponent.hh" // Implementation of class methods

#include "pylith/utils/JournalingComponent.hh" // USES JournalingComponent

#include "pylith/utils/journals.hh" // USES PYLITH_JOURNAL_*

// ----------------------------------------------------------------------
CPPUNIT_TEST_SUITE_REGISTRATION( pylith::utils::TestJournalingComponent );

// ----------------------------------------------------------------------
// Test constructor.
void
pylith::utils::TestJournalingComponent::testConstructor(void)
{ // testConstructor
    JournalingComponent journals;

    CPPUNIT_ASSERT_EQUAL(std::string(""), journals._name);
} // testConstructor

// ----------------------------------------------------------------------
// Test name().
void
pylith::utils::TestJournalingComponent::testName(void)
{ // testName
    JournalingComponent journals;
    CPPUNIT_ASSERT_EQUAL(std::string(""), std::string(journals.name()));

    const std::string& name = "my name";
    journals.name(name.c_str());
    CPPUNIT_ASSERT_EQUAL(name, std::string(journals.name()));
} // testName

namespace pylith {
    namespace utils {
        class TestJournals : public JournalingComponent {
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

// ----------------------------------------------------------------------
// Test PYLITH_JOURNAL_*.
void
pylith::utils::TestJournalingComponent::testJournals(void)
{ // testJournals
    const char* name = "test";
    journal::info_t info(name); info.activate();
    journal::debug_t debug(name); debug.activate();
    journal::warning_t warning(name); warning.activate();
    journal::error_t error(name); error.activate();

    TestJournals journals;
    journals.name("test");
    journals.test();
} // testJournals


// End of file
