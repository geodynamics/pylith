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

/**
 * @file unittests/libtests/utils/TestJournalingComponent.hh
 *
 * @brief C++ TestJournalingComponent object
 *
 * C++ unit testing for JournalingComponent.
 */

#if !defined(pylith_utils_testjournalingcomponent_hh)
#define pylith_utils_testjournalingcomponent_hh

#include <cppunit/extensions/HelperMacros.h>

/// Namespace for pylith package
namespace pylith {
    namespace utils {
        class TestJournalingComponent;
    } // utils
} // pylith

/// C++ unit testing for TestJournalingComponent
class pylith::utils::TestJournalingComponent : public CppUnit::TestFixture
{ // class TestJournalingComponent

// CPPUNIT TEST SUITE /////////////////////////////////////////////////
CPPUNIT_TEST_SUITE( TestJournalingComponent );

CPPUNIT_TEST( testConstructor );
CPPUNIT_TEST( testName );
CPPUNIT_TEST( testInitialize );
CPPUNIT_TEST( testAccessors );
CPPUNIT_TEST( testJournals );

CPPUNIT_TEST_SUITE_END();

// PUBLIC METHODS ///////////////////////////////////////////////////////
public:

/// Test constructor.
void testConstructor(void);

/// Test name().
void testName(void);

/// Test initialize().
void testInitialize(void);

/// Test accessors debug(), info(), error().
void testAccessors(void);

/// Test PYLITH_JOURNAL_DEBUG(), PYLITH_JOURNAL_INFO(), PYLITH_JOURNAL_ERROR().
void testJournals(void);

}; // class TestJournalingComponent

#endif // pylith_utils_testjournalingcomponent_hh


// End of file
