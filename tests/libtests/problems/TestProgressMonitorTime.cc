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

#include "TestProgressMonitorTime.hh" // Implementation of class methods

#include "pylith/problems/ProgressMonitorTime.hh" // USES ProgressMonitorTime

#include "pylith/utils/error.hh" // USES PYLITH_METHOD_BEGIN/END

// ---------------------------------------------------------------------------------------------------------------------
CPPUNIT_TEST_SUITE_REGISTRATION(pylith::problems::TestProgressMonitorTime);

// ---------------------------------------------------------------------------------------------------------------------
// Setup testing data.
void
pylith::problems::TestProgressMonitorTime::setUp(void) {
    _monitor = new ProgressMonitorTime();CPPUNIT_ASSERT(_monitor);
} // setUp


// ---------------------------------------------------------------------------------------------------------------------
// Tear down testing data.
void
pylith::problems::TestProgressMonitorTime::tearDown(void) {
    delete _monitor;_monitor = NULL;
} // tearDown


// ---------------------------------------------------------------------------------------------------------------------
// Test get/setUpdatePercent() and get/setFilename().
void
pylith::problems::TestProgressMonitorTime::testAccessors(void) {
    PYLITH_METHOD_BEGIN;

    CPPUNIT_ASSERT(_monitor);

    std::string unit = "second";
    CPPUNIT_ASSERT_EQUAL(unit, std::string(_monitor->getTimeUnit()));

    unit = "year";
    _monitor->setTimeUnit(unit.c_str());
    CPPUNIT_ASSERT_EQUAL(unit, std::string(_monitor->getTimeUnit()));

    PYLITH_METHOD_END;
} // testAccessors


// ---------------------------------------------------------------------------------------------------------------------
// Test open() and close().
void
pylith::problems::TestProgressMonitorTime::testOpenClose(void) {
    PYLITH_METHOD_BEGIN;

    CPPUNIT_ASSERT(_monitor);
    _monitor->open();
    _monitor->close();
    _monitor->close();

    PYLITH_METHOD_END;
} // testOpenClose


// ---------------------------------------------------------------------------------------------------------------------
// Test update().
void
pylith::problems::TestProgressMonitorTime::testUpdate(void) {
    PYLITH_METHOD_BEGIN;

    const double start = 1.0;
    const double stop = 11.0;
    const int numUpdateCalls = 7;
    const double current[numUpdateCalls] = { 1.0, 2.0, 2.1, 4.0, 6.1, 6.2, 10.0 };
    const int numUpdates = 5;
    const double percentComplete[numUpdates] = { 0.0, 10.0, 30.0, 51.0, 90.0 };
    const double tolerance = 1.0e-6;
    const char* filename = "progress_time.txt";

    CPPUNIT_ASSERT(_monitor);
    _monitor->setFilename(filename);
    _monitor->open();
    for (int i = 0; i < numUpdateCalls; ++i) {
        _monitor->update(current[i], start, stop);
    } // for
    _monitor->close();

    // Check output
    std::ifstream fin(filename);
    CPPUNIT_ASSERT_MESSAGE("Could not open progress monitor file.", fin.is_open());

    int count = 0;
    const int maxlen = 1024;
    char buffer[maxlen];
    fin.getline(buffer, maxlen); // Ignore header
    fin.getline(buffer, maxlen);
    while (fin.good()) {
        const double percentCompleteValue = stof(std::string(buffer).substr(43, 12));
        CPPUNIT_ASSERT_DOUBLES_EQUAL_MESSAGE("Incorrect percent complete in progress monitor output.",
                                             percentComplete[count++], percentCompleteValue, tolerance);
        fin.getline(buffer, maxlen);
    } // while
    CPPUNIT_ASSERT_EQUAL_MESSAGE("Incorrect number of updates in progress monitor output.", count, numUpdates);
    fin.close();

    PYLITH_METHOD_END;
} // testUpdate


// End of file
