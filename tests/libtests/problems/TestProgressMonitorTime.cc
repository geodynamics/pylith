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

#include "pylith/utils/GenericComponent.hh" // ISA GenericComponent

#include "pylith/problems/ProgressMonitorTime.hh" // USES ProgressMonitorTime

#include "pylith/utils/error.hh" // USES PYLITH_METHOD_BEGIN/END

#include "catch2/catch_test_macros.hpp"
#include "catch2/matchers/catch_matchers_floating_point.hpp"

/// Namespace for pylith package
namespace pylith {
    namespace problems {
        class TestProgressMonitorTime;
    } // problems
} // pylith

// ------------------------------------------------------------------------------------------------
class pylith::problems::TestProgressMonitorTime : public pylith::utils::GenericComponent {
public:

    /// Constructor
    TestProgressMonitorTime(void);

    /// Destructor.
    ~TestProgressMonitorTime(void);

    /// Test get/setTimeUnit().
    void testAccessors(void);

    /// Test open() and close().
    void testOpenClose(void);

    /// Test update().
    void testUpdate(void);

private:

    pylith::problems::ProgressMonitorTime* _monitor; ///< Test subject.

}; // class TestProgressMonitorTime

// ------------------------------------------------------------------------------------------------
TEST_CASE("TestProgressMonitorTime::testAccessors", "[TestProgressMonitorTime]") {
    pylith::problems::TestProgressMonitorTime().testAccessors();
}
TEST_CASE("TestProgressMonitorTime::testOpenClose", "[TestProgressMonitorTime]") {
    pylith::problems::TestProgressMonitorTime().testOpenClose();
}
TEST_CASE("TestProgressMonitorTime::testUpdate", "[TestProgressMonitorTime]") {
    pylith::problems::TestProgressMonitorTime().testUpdate();
}

// ------------------------------------------------------------------------------------------------
// Setup testing data.
pylith::problems::TestProgressMonitorTime::TestProgressMonitorTime(void) {
    _monitor = new ProgressMonitorTime();assert(_monitor);
} // setUp


// ------------------------------------------------------------------------------------------------
// Tear down testing data.
pylith::problems::TestProgressMonitorTime::~TestProgressMonitorTime(void) {
    delete _monitor;_monitor = NULL;
} // tearDown


// ------------------------------------------------------------------------------------------------
// Test get/setUpdatePercent() and get/setFilename().
void
pylith::problems::TestProgressMonitorTime::testAccessors(void) {
    PYLITH_METHOD_BEGIN;
    assert(_monitor);

    std::string unit = "second";
    CHECK(unit == std::string(_monitor->getTimeUnit()));

    unit = "year";
    _monitor->setTimeUnit(unit.c_str());
    CHECK(unit == std::string(_monitor->getTimeUnit()));

    PYLITH_METHOD_END;
} // testAccessors


// ------------------------------------------------------------------------------------------------
// Test open() and close().
void
pylith::problems::TestProgressMonitorTime::testOpenClose(void) {
    PYLITH_METHOD_BEGIN;
    assert(_monitor);

    _monitor->open();
    _monitor->close();
    _monitor->close();

    PYLITH_METHOD_END;
} // testOpenClose


// ------------------------------------------------------------------------------------------------
// Test update().
void
pylith::problems::TestProgressMonitorTime::testUpdate(void) {
    PYLITH_METHOD_BEGIN;
    assert(_monitor);

    const double start = 1.0;
    const double stop = 11.0;
    const int numUpdateCalls = 7;
    const double current[numUpdateCalls] = { 1.0, 2.0, 2.1, 4.0, 6.1, 6.2, 10.0 };
    const int numUpdates = 5;
    const double percentComplete[numUpdates] = { 0.0, 10.0, 30.0, 51.0, 90.0 };
    const double tolerance = 1.0e-6;
    const char* filename = "progress_time.txt";

    _monitor->setFilename(filename);
    _monitor->open();
    for (int i = 0; i < numUpdateCalls; ++i) {
        _monitor->update(current[i], start, stop);
    } // for
    _monitor->close();

    // Check output
    std::ifstream fin(filename);
    REQUIRE(fin.is_open());

    int count = 0;
    const int maxlen = 1024;
    char buffer[maxlen];
    fin.getline(buffer, maxlen); // Ignore header
    fin.getline(buffer, maxlen);
    while (fin.good()) {
        const double percentCompleteValue = stof(std::string(buffer).substr(43, 12));
        CHECK_THAT(percentCompleteValue, Catch::Matchers::WithinAbs(percentComplete[count++], tolerance));
        fin.getline(buffer, maxlen);
    } // while
    CHECK(count == numUpdates);
    fin.close();

    PYLITH_METHOD_END;
} // testUpdate


// End of file
