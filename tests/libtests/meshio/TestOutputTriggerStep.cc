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

#include "TestOutputTriggerStep.hh" // Implementation of class methods

#include "pylith/meshio/OutputTriggerStep.hh" // USES OutputTriggerStep

// ---------------------------------------------------------------------------------------------------------------------
CPPUNIT_TEST_SUITE_REGISTRATION(pylith::meshio::TestOutputTriggerStep);

// ---------------------------------------------------------------------------------------------------------------------
// Test setNumStepsSkip() and getNumStepsSkip().
void
pylith::meshio::TestOutputTriggerStep::testNumStepsSkip(void) {
    OutputTriggerStep trigger;

    int numSkip = 0; // default
    CPPUNIT_ASSERT_EQUAL(numSkip, trigger.getNumStepsSkip());

    numSkip = 2;
    trigger.setNumStepsSkip(numSkip);
    CPPUNIT_ASSERT_EQUAL(numSkip, trigger.getNumStepsSkip());
} // testNumStepsSkip


// ---------------------------------------------------------------------------------------------------------------------
// Test shouldWrite().
void
pylith::meshio::TestOutputTriggerStep::testShouldWrite(void) {
    OutputTriggerStep trigger;

    const PylithReal dt = 0.1;
    PylithReal t = 0.0;
    PylithInt tindex = 0;
    CPPUNIT_ASSERT(trigger.shouldWrite(t, tindex++));t += dt;
    CPPUNIT_ASSERT(trigger.shouldWrite(t, tindex++));t += dt;
    CPPUNIT_ASSERT(trigger.shouldWrite(t, tindex++));t += dt;

    trigger.setNumStepsSkip(1);
    CPPUNIT_ASSERT(!trigger.shouldWrite(t, tindex++));t += dt;
    CPPUNIT_ASSERT(trigger.shouldWrite(t, tindex++));t += dt;
    CPPUNIT_ASSERT(!trigger.shouldWrite(t, tindex++));t += dt;
    CPPUNIT_ASSERT(trigger.shouldWrite(t, tindex++));t += dt;
    CPPUNIT_ASSERT(!trigger.shouldWrite(t, tindex++));t += dt;

    trigger.setNumStepsSkip(2);
    CPPUNIT_ASSERT(!trigger.shouldWrite(t, tindex++));t += dt;
    CPPUNIT_ASSERT(trigger.shouldWrite(t, tindex++));t += dt;
    CPPUNIT_ASSERT(!trigger.shouldWrite(t, tindex++));t += dt;
    CPPUNIT_ASSERT(!trigger.shouldWrite(t, tindex++));t += dt;
    CPPUNIT_ASSERT(trigger.shouldWrite(t, tindex++));t += dt;
    CPPUNIT_ASSERT(!trigger.shouldWrite(t, tindex++));t += dt;
    CPPUNIT_ASSERT(!trigger.shouldWrite(t, tindex++));t += dt;
    CPPUNIT_ASSERT(trigger.shouldWrite(t, tindex++));t += dt;
    CPPUNIT_ASSERT(!trigger.shouldWrite(t, tindex++));t += dt;
} // testShouldWrite


// End of file
