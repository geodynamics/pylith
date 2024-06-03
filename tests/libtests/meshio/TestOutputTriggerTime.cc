// =================================================================================================
// This code is part of PyLith, developed through the Computational Infrastructure
// for Geodynamics (https://github.com/geodynamics/pylith).
//
// Copyright (c) 2010-2024, University of California, Davis and the PyLith Development Team.
// All rights reserved.
//
// See https://mit-license.org/ and LICENSE.md and for license information.
// =================================================================================================

#include <portinfo>

#include "pylith/utils/GenericComponent.hh" // ISA GenericComponent

#include "pylith/meshio/OutputTriggerTime.hh" // USES OutputTriggerTime

#include "catch2/catch_test_macros.hpp"
#include "catch2/matchers/catch_matchers_floating_point.hpp"

// ------------------------------------------------------------------------------------------------
namespace pylith {
    namespace meshio {
        class TestOutputTriggerTime;
    } // meshio
} // pylith

// ------------------------------------------------------------------------------------------------
class pylith::meshio::TestOutputTriggerTime : public pylith::utils::GenericComponent {
    // PUBLIC METHODS /////////////////////////////////////////////////////////////////////////////
public:

    /// Test setTimeSkip() and getTimeSkip().
    static
    void testTimeSkip(void);

    /// Test shouldWrite().
    static
    void testShouldWrite(void);

}; // TestOutputTriggerTime

// ------------------------------------------------------------------------------------------------
TEST_CASE("TestOutputTriggerTime::testTimeSkip", "[TestOutputTriggerTime][testTimeSkip]") {
    pylith::meshio::TestOutputTriggerTime::testTimeSkip();
}
TEST_CASE("TestOutputTriggerTime::testShouldWrite", "[TestOutputTriggerTime][testShouldWrite]") {
    pylith::meshio::TestOutputTriggerTime::testShouldWrite();
}

// ------------------------------------------------------------------------------------------------
// Test setTimeSkip() and getTimeSkip().
void
pylith::meshio::TestOutputTriggerTime::testTimeSkip(void) {
    const double tolerance = 1.0e-6;

    OutputTriggerTime trigger;

    PylithReal tskip = 0.0; // default
    CHECK_THAT(trigger.getTimeSkip(), Catch::Matchers::WithinAbs(tskip, tolerance));

    tskip = 1.5;
    trigger.setTimeSkip(tskip);
    CHECK_THAT(trigger.getTimeSkip(), Catch::Matchers::WithinAbs(tskip, tolerance));
} // testTimeSkip


// ------------------------------------------------------------------------------------------------
// Test shouldWrite().
void
pylith::meshio::TestOutputTriggerTime::testShouldWrite(void) {
    OutputTriggerTime trigger;

    const PylithReal dt = 0.1;
    PylithReal t = 0.0;
    PylithInt tindex = 0;
    CHECK(true == trigger.shouldWrite(t, tindex++));t += dt;
    CHECK(true == trigger.shouldWrite(t, tindex++));t += dt;
    CHECK(true == trigger.shouldWrite(t, tindex++));t += dt;

    trigger.setTimeSkip(0.1999);
    CHECK(false == trigger.shouldWrite(t, tindex++));t += dt;
    CHECK(true == trigger.shouldWrite(t, tindex++));t += dt;
    CHECK(false == trigger.shouldWrite(t, tindex++));t += dt;
    CHECK(true == trigger.shouldWrite(t, tindex++));t += dt;
    CHECK(false == trigger.shouldWrite(t, tindex++));t += dt;

    trigger.setTimeSkip(0.2999);
    CHECK(false == trigger.shouldWrite(t, tindex++));t += dt;
    CHECK(true == trigger.shouldWrite(t, tindex++));t += dt;
    CHECK(false == trigger.shouldWrite(t, tindex++));t += dt;
    CHECK(false == trigger.shouldWrite(t, tindex++));t += dt;
    CHECK(true == trigger.shouldWrite(t, tindex++));t += dt;
    CHECK(false == trigger.shouldWrite(t, tindex++));t += dt;
    CHECK(false == trigger.shouldWrite(t, tindex++));t += dt;
    CHECK(true == trigger.shouldWrite(t, tindex++));t += dt;
    CHECK(false == trigger.shouldWrite(t, tindex++));t += dt;
} // testShouldWrite


// End of file
