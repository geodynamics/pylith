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
    CPPUNIT_ASSERT(!journals._debug);
    CPPUNIT_ASSERT(!journals._info);
    CPPUNIT_ASSERT(!journals._error);
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

// ----------------------------------------------------------------------
// Test initialize().
void
pylith::utils::TestJournalingComponent::testInitialize(void)
{ // testInitialize
    JournalingComponent journals;
    journals.name("my class");

    journals.initialize();

    CPPUNIT_ASSERT(journals._debug);
    CPPUNIT_ASSERT(journals._info);
    CPPUNIT_ASSERT(journals._error);
} // testInitialize

// ----------------------------------------------------------------------
// Test debug(), info(), error().
void
pylith::utils::TestJournalingComponent::testAccessors(void)
{ // testAccessors
    JournalingComponent journals;
    journals.name("TestJournalingComponent");
    journals.initialize();

    journal::debug_t& debug = journals.debug();
    debug << journal::at(__HERE__) << "CORRECT: This is a debug message." << journal::endl;

    journal::info_t& info = journals.info();
    info << journal::at(__HERE__) << "CORRECT: This is a info message." << journal::endl;

    journal::error_t& error = journals.error();
    error << journal::at(__HERE__) << "CORRECT: This is a error message." << journal::endl;

} // testAccessors

namespace pylith {
    namespace utils {
        class TestJournals : public JournalingComponent {
public:
        void test(void) {
            debug().activate();
            PYLITH_JOURNAL_DEBUG("CORRECT: This is a debug message.");

            info().activate();
            PYLITH_JOURNAL_INFO("CORRECT: This is an info mesasge.");

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
    TestJournals journals;
    journals.name("test");
    journals.initialize();
    journals.test();
} // testJournals


// End of file
