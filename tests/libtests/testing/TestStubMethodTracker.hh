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
 * @file tests/libtests/utils/TestStubMethodTracker.hh
 *
 * @brief C++ TestStubMethodTracker object
 *
 * C++ unit testing for StubMethodTracker.
 */

#if !defined(pylith_utils_teststubmethodtracker_hh)
#define pylith_utils_teststubmethodtracker_hh

#include <cppunit/extensions/HelperMacros.h>

/// Namespace for pylith package
namespace pylith {
    namespace testing {
        class TestStubMethodTracker;
    } // utils
} // pylith

/// C++ unit testing for TestStubMethodTracker
class pylith::testing::TestStubMethodTracker : public CppUnit::TestFixture {
    // CPPUNIT TEST SUITE /////////////////////////////////////////////////
    CPPUNIT_TEST_SUITE(TestStubMethodTracker);

    CPPUNIT_TEST(testMethodCount);

    CPPUNIT_TEST_SUITE_END();

    // PUBLIC METHODS ///////////////////////////////////////////////////////
public:

    /// Test methodCount().
    void testMethodCount(void);

}; // class TestStubMethodTracker

#endif // pylith_utils_teststubmethodtracker_hh

// End of file
