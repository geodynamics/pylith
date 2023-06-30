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

#include "tests/src/ProgressMonitorStub.hh" // USES ProgressMonitorStub
#include "tests/src/StubMethodTracker.hh" // USES StubMethodTracker

#include "pylith/utils/error.hh" // USES PYLITH_METHOD_BEGIN/END

#include "catch2/catch_test_macros.hpp"
#include "catch2/matchers/catch_matchers_floating_point.hpp"

// ------------------------------------------------------------------------------------------------
namespace pylith {
    namespace problems {
        class TestProgressMonitor;
    } // problems
} // pylith

class pylith::problems::TestProgressMonitor : public pylith::utils::GenericComponent {
public:

    /// Constructor.
    TestProgressMonitor(void);

    /// Destructor.
    ~TestProgressMonitor(void);

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
TEST_CASE("TestProgressMonitor::testAccessors", "[TestProgressMonitor]") {
    pylith::problems::TestProgressMonitor().testAccessors();
}
TEST_CASE("TestProgressMonitor::testOpenClose", "[TestProgressMonitor]") {
    pylith::problems::TestProgressMonitor().testOpenClose();
}
TEST_CASE("TestProgressMonitor::testUpdate", "[TestProgressMonitor]") {
    pylith::problems::TestProgressMonitor().testUpdate();
}

// ------------------------------------------------------------------------------------------------
// Constructor.
pylith::problems::TestProgressMonitor::TestProgressMonitor(void) {
    _monitor = new ProgressMonitorStub();assert(_monitor);
} // setUp


// ------------------------------------------------------------------------------------------------
// Destructor.
pylith::problems::TestProgressMonitor::~TestProgressMonitor(void) {
    delete _monitor;_monitor = NULL;
} // tearDown


// ------------------------------------------------------------------------------------------------
// Test get/setUpdatePercent() and get/setFilename().
void
pylith::problems::TestProgressMonitor::testAccessors(void) {
    PYLITH_METHOD_BEGIN;
    assert(_monitor);
    const double tolerance = 1.0e-6;

    { // updatePercent
        const double valueDefault = 5.0;
        CHECK_THAT(_monitor->getUpdatePercent(), Catch::Matchers::WithinAbs(valueDefault, tolerance));

        const double value = 2.0;
        _monitor->setUpdatePercent(value);
        CHECK_THAT(_monitor->getUpdatePercent(), Catch::Matchers::WithinAbs(value, tolerance));
    } // updatePercent

    { // filename
        const std::string& valueDefault = "progress.txt";
        CHECK(valueDefault == std::string(_monitor->getFilename()));

        const std::string& value = "progress2.txt";
        _monitor->setFilename(value.c_str());
        CHECK(value == std::string(_monitor->getFilename()));
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

    assert(_monitor);
    _monitor->open();
    CHECK(size_t(1) == tracker.getMethodCount("pylith::problems::ProgressMonitorStub::_open"));

    _monitor->close();
    CHECK(size_t(1) == tracker.getMethodCount("pylith::problems::ProgressMonitorStub::_close"));

    _monitor->close();
    CHECK(size_t(2) == tracker.getMethodCount("pylith::problems::ProgressMonitorStub::_close"));

    PYLITH_METHOD_END;
} // testOpenClose


// ------------------------------------------------------------------------------------------------
// Test update().
void
pylith::problems::TestProgressMonitor::testUpdate(void) {
    PYLITH_METHOD_BEGIN;

    pylith::testing::StubMethodTracker tracker;
    tracker.clear();

    assert(_monitor);
    _monitor->open();
    _monitor->close();

    PYLITH_METHOD_END;
} // testUpdate


// End of file
