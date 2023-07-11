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

#include "pylith/utils/GenericComponent.hh" // ISA GenericComponent

#include "tests/src/StubMethodTracker.hh" // USES StubMethodTracker

#include "catch2/catch_test_macros.hpp"

// ------------------------------------------------------------------------------------------------
namespace pylith {
    namespace testing {
        class TestStubMethodTracker;
    }
}

class pylith::testing::TestStubMethodTracker : public pylith::utils::GenericComponent {
public:

    /// Test methodCount().
    static
    void testMethodCount(void);

}; // class TestStubMethodTracker

// ------------------------------------------------------------------------------------------------
TEST_CASE("TestStubMethodTracker::testMethodCount", "[TestStubMethodTracker]") {
    pylith::testing::TestStubMethodTracker::testMethodCount();
}

// ------------------------------------------------------------------------------------------------
void
pylith::testing::TestStubMethodTracker::testMethodCount(void) {
    StubMethodTracker tracker;

    const char* methodA = "pylith::testing::A";
    const char* methodB = "pylith::testing::B";
    CHECK(size_t(0) == tracker.getMethodCount(methodA));
    CHECK(size_t(0) == tracker.getMethodCount(methodA));
    CHECK(size_t(0) == tracker.getMethodCount(methodB));

    tracker.methodCalled(methodA);
    CHECK(size_t(1) == tracker.getMethodCount(methodA));
    CHECK(size_t(0) == tracker.getMethodCount(methodB));

    tracker.methodCalled(methodB);
    CHECK(size_t(1) == tracker.getMethodCount(methodA));
    CHECK(size_t(1) == tracker.getMethodCount(methodB));

    tracker.methodCalled(methodA);
    CHECK(size_t(2) == tracker.getMethodCount(methodA));
    CHECK(size_t(1) == tracker.getMethodCount(methodB));

    StubMethodTracker trackerB(methodB);
    CHECK(size_t(2) == tracker.getMethodCount(methodA));
    CHECK(size_t(2) == tracker.getMethodCount(methodB));

    tracker.clear();
    CHECK(size_t(0) == tracker.getMethodCount(methodA));
    CHECK(size_t(0) == tracker.getMethodCount(methodB));
} // testMethodCount


// End of file
