// =================================================================================================
// This code is part of PyLith, developed through the Computational Infrastructure
// for Geodynamics (https://github.com/geodynamics/pylith).
//
// Copyright (c) 2010-2023, University of California, Davis and the PyLith Development Team.
// All rights reserved.
//
// See https://mit-license.org/ and LICENSE.md and for license information.
// =================================================================================================

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
