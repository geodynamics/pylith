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
/** C++ unit testing for ProgressMonitor.
 */

#include <portinfo>

#include <cppunit/extensions/HelperMacros.h>

#include "pylith/testing/ProgressMonitorStub.hh" // USES ProgressMonitorStub
#include "pylith/testing/StubMethodTracker.hh" // USES StubMethodTracker

#include "pylith/utils/error.hh" // USES PYLITH_METHOD_BEGIN/END

// ------------------------------------------------------------------------------------------------
/// Namespace for pylith package
namespace pylith {
    namespace problems {
        class TestProgressMonitor;
    } // problems
} // pylith

class pylith::problems::TestProgressMonitor : public CppUnit::TestFixture {
    CPPUNIT_TEST_SUITE(TestProgressMonitor);

    CPPUNIT_TEST(testAccessors);
    CPPUNIT_TEST(testOpenClose);
    CPPUNIT_TEST(testUpdate);

    CPPUNIT_TEST_SUITE_END();

public:

    /// Setup testing data.
    void setUp(void);

    /// Tear down testing data.
    void tearDown(void);

    /// Test get/setUpdatePercent() and get/setFilename().
    void testAccessors(void);

    /// Test open() and close().
    void testOpenClose(void);

    /// Test update().
    void testUpdate(void);

private:

    pylith::problems::ProgressMonitorStub* _monitor; ///< Test subject.

}; // class TestProgressMonitor

// ------------------------------------------------------------------------------------------------
CPPUNIT_TEST_SUITE_REGISTRATION(pylith::problems::TestProgressMonitor);

// ------------------------------------------------------------------------------------------------
// Setup testing data.
void
pylith::problems::TestProgressMonitor::setUp(void) {
    _monitor = new ProgressMonitorStub();CPPUNIT_ASSERT(_monitor);
} // setUp


// ------------------------------------------------------------------------------------------------
// Tear down testing data.
void
pylith::problems::TestProgressMonitor::tearDown(void) {
    delete _monitor;_monitor = NULL;
} // tearDown


// ------------------------------------------------------------------------------------------------
// Test get/setUpdatePercent() and get/setFilename().
void
pylith::problems::TestProgressMonitor::testAccessors(void) {
    PYLITH_METHOD_BEGIN;

    CPPUNIT_ASSERT(_monitor);
    const double tolerance = 1.0e-6;

    { // updatePercent
        const double valueDefault = 5.0;
        CPPUNIT_ASSERT_DOUBLES_EQUAL(valueDefault, _monitor->getUpdatePercent(), tolerance);

        const double value = 2.0;
        _monitor->setUpdatePercent(value);
        CPPUNIT_ASSERT_DOUBLES_EQUAL(value, _monitor->getUpdatePercent(), tolerance);
    } // updatePercent

    { // filename
        const std::string& valueDefault = "progress.txt";
        CPPUNIT_ASSERT_EQUAL(valueDefault, std::string(_monitor->getFilename()));

        const std::string& value = "progress2.txt";
        _monitor->setFilename(value.c_str());
        CPPUNIT_ASSERT_EQUAL(value, std::string(_monitor->getFilename()));
    } // filename

    PYLITH_METHOD_END;
} // testAccessors


// ------------------------------------------------------------------------------------------------
// Test open() and close().
void
pylith::problems::TestProgressMonitor::testOpenClose(void) {
    PYLITH_METHOD_BEGIN;

    pylith::testing::StubMethodTracker tracker;
    tracker.clear();

    CPPUNIT_ASSERT(_monitor);
    _monitor->open();
    CPPUNIT_ASSERT_EQUAL(size_t(1), tracker.getMethodCount("pylith::problems::ProgressMonitorStub::_open"));

    _monitor->close();
    CPPUNIT_ASSERT_EQUAL(size_t(1), tracker.getMethodCount("pylith::problems::ProgressMonitorStub::_close"));

    _monitor->close();
    CPPUNIT_ASSERT_EQUAL(size_t(2), tracker.getMethodCount("pylith::problems::ProgressMonitorStub::_close"));

    PYLITH_METHOD_END;
} // testOpenClose


// ------------------------------------------------------------------------------------------------
// Test update().
void
pylith::problems::TestProgressMonitor::testUpdate(void) {
    PYLITH_METHOD_BEGIN;

    pylith::testing::StubMethodTracker tracker;
    tracker.clear();

    CPPUNIT_ASSERT(_monitor);
    _monitor->open();
    _monitor->close();

    PYLITH_METHOD_END;
} // testUpdate


// End of file
