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

#include "TestProgressMonitor.hh" // Implementation of class methods

#include "pylith/testing/ProgressMonitorStub.hh" // USES ProgressMonitorStub
#include "pylith/testing/StubMethodTracker.hh" // USES StubMethodTracker

#include "pylith/utils/error.hh" // USES PYLITH_METHOD_BEGIN/END

// ---------------------------------------------------------------------------------------------------------------------
CPPUNIT_TEST_SUITE_REGISTRATION(pylith::problems::TestProgressMonitor);

// ---------------------------------------------------------------------------------------------------------------------
// Setup testing data.
void
pylith::problems::TestProgressMonitor::setUp(void) {
    _monitor = new ProgressMonitorStub();CPPUNIT_ASSERT(_monitor);
} // setUp


// ---------------------------------------------------------------------------------------------------------------------
// Tear down testing data.
void
pylith::problems::TestProgressMonitor::tearDown(void) {
    delete _monitor;_monitor = NULL;
} // tearDown


// ---------------------------------------------------------------------------------------------------------------------
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


// ---------------------------------------------------------------------------------------------------------------------
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


// ---------------------------------------------------------------------------------------------------------------------
// Test update().
void
pylith::problems::TestProgressMonitor::testUpdate(void) {
    PYLITH_METHOD_BEGIN;

    pylith::testing::StubMethodTracker tracker;
    tracker.clear();

    CPPUNIT_ASSERT(_monitor);
    _monitor->open();

    size_t count = 0;
    double current = 2.0;
    double start = 2.0;
    double stop = 12.0;
    const double tolerance = 1.0e-6;

    _monitor->update(current, start, stop);
    CPPUNIT_ASSERT_EQUAL(++count, tracker.getMethodCount("pylith::problems::ProgressMonitorStub::_update"));
    CPPUNIT_ASSERT_DOUBLES_EQUAL(current, _monitor->_state.current, tolerance);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(0.0, _monitor->_state.percentComplete, tolerance);

    current = 4.0;
    _monitor->update(current, start, stop);
    CPPUNIT_ASSERT_EQUAL(++count, tracker.getMethodCount("pylith::problems::ProgressMonitorStub::_update"));
    CPPUNIT_ASSERT_DOUBLES_EQUAL(current, _monitor->_state.current, tolerance);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(20.0, _monitor->_state.percentComplete, tolerance);

    current = 4.1;
    _monitor->update(current, start, stop);
    CPPUNIT_ASSERT_EQUAL(count, tracker.getMethodCount("pylith::problems::ProgressMonitorStub::_update"));
    CPPUNIT_ASSERT_DOUBLES_EQUAL(4.0, _monitor->_state.current, tolerance);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(20.0, _monitor->_state.percentComplete, tolerance);

    _monitor->close();

    PYLITH_METHOD_END;
} // testUpdate


// End of file
