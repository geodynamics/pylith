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

#include "pylith/problems/ProgressMonitorStep.hh" // USES ProgressMonitorStep

#include "pylith/utils/error.hh" // USES PYLITH_METHOD_BEGIN/END

#include "catch2/catch_test_macros.hpp"
#include "catch2/matchers/catch_matchers_floating_point.hpp"

// ------------------------------------------------------------------------------------------------
/// Namespace for pylith package
namespace pylith {
    namespace problems {
        class TestProgressMonitorStep;
    } // problems
} // pylith

class pylith::problems::TestProgressMonitorStep : public pylith::utils::GenericComponent {
public:

    /// Constructor.
    TestProgressMonitorStep(void);

    /// Destructor.
    ~TestProgressMonitorStep(void);

    /// Test open() and close().
    void testOpenClose(void);

    /// Test update().
    void testUpdate(void);

private:

    pylith::problems::ProgressMonitorStep* _monitor; ///< Test subject.

}; // class TestProgressMonitorStep

// ------------------------------------------------------------------------------------------------
TEST_CASE("TestProgressMonitorStep::testOpenClose", "[TestProgressMonitorStep]") {
    pylith::problems::TestProgressMonitorStep().testOpenClose();
}
TEST_CASE("TestProgressMonitorStep::testUpdate", "[TestProgressMonitorStep]") {
    pylith::problems::TestProgressMonitorStep().testUpdate();
}

// ------------------------------------------------------------------------------------------------
// Constructor.
pylith::problems::TestProgressMonitorStep::TestProgressMonitorStep(void) {
    _monitor = new ProgressMonitorStep();assert(_monitor);
} // setUp


// ------------------------------------------------------------------------------------------------
// Destructor.
pylith::problems::TestProgressMonitorStep::~TestProgressMonitorStep(void) {
    delete _monitor;_monitor = NULL;
} // tearDown


// ------------------------------------------------------------------------------------------------
// Test open() and close().
void
pylith::problems::TestProgressMonitorStep::testOpenClose(void) {
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
pylith::problems::TestProgressMonitorStep::testUpdate(void) {
    PYLITH_METHOD_BEGIN;
    assert(_monitor);

    const size_t start = 2;
    const size_t stop = 52;
    const int numUpdateCalls = 7;
    const size_t current[numUpdateCalls] = { 2, 5, 7, 8, 10, 11, 15 };
    const int numUpdates = 5;
    const double percentComplete[numUpdates] = { 0.0, 6.0, 10.0, 16.0, 26.0 };
    const double tolerance = 1.0e-6;
    const char* filename = "progress_step.txt";

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
        const double percentCompleteValue = stof(std::string(buffer).substr(41, 12));
        CHECK_THAT(percentCompleteValue, Catch::Matchers::WithinAbs(percentComplete[count++], tolerance));
        fin.getline(buffer, maxlen);
    } // while
    CHECK(count == numUpdates);
    fin.close();

    PYLITH_METHOD_END;
} // testUpdate


// End of file
