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

/**
 * @file tests/libtests/utils/TestPyreComponent.hh
 *
 * @brief C++ TestPyreComponent object
 *
 * C++ unit testing for PyreComponent.
 */

#if !defined(pylith_utils_testpyrecomponent_hh)
#define pylith_utils_testpyrecomponent_hh

#include <cppunit/extensions/HelperMacros.h>

/// Namespace for pylith package
namespace pylith {
    namespace utils {
        class TestPyreComponent;
    } // utils
} // pylith

/// C++ unit testing for TestPyreComponent
class pylith::utils::TestPyreComponent : public CppUnit::TestFixture {

    // CPPUNIT TEST SUITE /////////////////////////////////////////////////
    CPPUNIT_TEST_SUITE(TestPyreComponent);

    CPPUNIT_TEST(testConstructor);
    CPPUNIT_TEST(testName);
    CPPUNIT_TEST(testIdentifier);
    CPPUNIT_TEST(testJournals);

    CPPUNIT_TEST_SUITE_END();

    // PUBLIC METHODS ///////////////////////////////////////////////////////
public:

    /// Test constructor.
    void testConstructor(void);

    /// Test name().
    void testName(void);

    /// Test identifier().
    void testIdentifier(void);

    /// Test PYLITH_JOURNAL_DEBUG(), PYLITH_JOURNAL_INFO(), PYLITH_JOURNAL_ERROR().
    void testJournals(void);

}; // class TestPyreComponent

#endif // pylith_utils_testpyrecomponent_hh


// End of file
