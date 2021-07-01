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

#include "TestStubMethodTracker.hh" // Implementation of class methods

#include "pylith/testing/StubMethodTracker.hh" // USES StubMethodTracker

// ----------------------------------------------------------------------
CPPUNIT_TEST_SUITE_REGISTRATION(pylith::testing::TestStubMethodTracker);

// ----------------------------------------------------------------------
// Test constructor.
void
pylith::testing::TestStubMethodTracker::testMethodCount(void) {
    StubMethodTracker tracker;

    const char* methodA = "pylith::testing::A";
    const char* methodB = "pylith::testing::B";
    CPPUNIT_ASSERT_EQUAL(size_t(0), tracker.getMethodCount(methodA));
    CPPUNIT_ASSERT_EQUAL(size_t(0), tracker.getMethodCount(methodA));
    CPPUNIT_ASSERT_EQUAL(size_t(0), tracker.getMethodCount(methodB));

    tracker.methodCalled(methodA);
    CPPUNIT_ASSERT_EQUAL(size_t(1), tracker.getMethodCount(methodA));
    CPPUNIT_ASSERT_EQUAL(size_t(0), tracker.getMethodCount(methodB));

    tracker.methodCalled(methodB);
    CPPUNIT_ASSERT_EQUAL(size_t(1), tracker.getMethodCount(methodA));
    CPPUNIT_ASSERT_EQUAL(size_t(1), tracker.getMethodCount(methodB));

    tracker.methodCalled(methodA);
    CPPUNIT_ASSERT_EQUAL(size_t(2), tracker.getMethodCount(methodA));
    CPPUNIT_ASSERT_EQUAL(size_t(1), tracker.getMethodCount(methodB));

    StubMethodTracker trackerB(methodB);
    CPPUNIT_ASSERT_EQUAL(size_t(2), tracker.getMethodCount(methodA));
    CPPUNIT_ASSERT_EQUAL(size_t(2), tracker.getMethodCount(methodB));

    tracker.clear();
    CPPUNIT_ASSERT_EQUAL(size_t(0), tracker.getMethodCount(methodA));
    CPPUNIT_ASSERT_EQUAL(size_t(0), tracker.getMethodCount(methodB));
} // testMethodCount


// End of file
