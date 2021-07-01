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
 * @file tests/libtests/utils/TestGenericComponent.hh
 *
 * @brief C++ TestGenericComponent object
 *
 * C++ unit testing for GenericComponent.
 */

#if !defined(pylith_utils_testgenericcomponent_hh)
#define pylith_utils_testgenericcomponent_hh

#include <cppunit/extensions/HelperMacros.h>

/// Namespace for pylith package
namespace pylith {
    namespace utils {
        class TestGenericComponent;
    } // utils
} // pylith

/// C++ unit testing for TestGenericComponent
class pylith::utils::TestGenericComponent : public CppUnit::TestFixture {

    // CPPUNIT TEST SUITE /////////////////////////////////////////////////
    CPPUNIT_TEST_SUITE(TestGenericComponent);

    CPPUNIT_TEST(testConstructor);
    CPPUNIT_TEST(testName);
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

}; // class TestGenericComponent

#endif // pylith_utils_testgenericcomponent_hh


// End of file
