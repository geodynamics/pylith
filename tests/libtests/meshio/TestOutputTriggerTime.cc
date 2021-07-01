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

#include "TestOutputTriggerTime.hh" // Implementation of class methods

#include "pylith/meshio/OutputTriggerTime.hh" // USES OutputTriggerTime

// ---------------------------------------------------------------------------------------------------------------------
CPPUNIT_TEST_SUITE_REGISTRATION(pylith::meshio::TestOutputTriggerTime);

// ---------------------------------------------------------------------------------------------------------------------
// Test setTimeSkip() and getTimeSkip().
void
pylith::meshio::TestOutputTriggerTime::testTimeSkip(void) {
    OutputTriggerTime trigger;

    PylithReal tskip = 0.0; // default
    CPPUNIT_ASSERT_EQUAL(tskip, trigger.getTimeSkip());

    tskip = 1.5;
    trigger.setTimeSkip(tskip);
    CPPUNIT_ASSERT_EQUAL(tskip, trigger.getTimeSkip());
} // testTimeSkip


// ---------------------------------------------------------------------------------------------------------------------
// Test shouldWrite().
void
pylith::meshio::TestOutputTriggerTime::testShouldWrite(void) {
    OutputTriggerTime trigger;

    const PylithReal dt = 0.1;
    PylithReal t = 0.0;
    PylithInt tindex = 0;
    CPPUNIT_ASSERT(trigger.shouldWrite(t, tindex++));t += dt;
    CPPUNIT_ASSERT(trigger.shouldWrite(t, tindex++));t += dt;
    CPPUNIT_ASSERT(trigger.shouldWrite(t, tindex++));t += dt;

    trigger.setTimeSkip(0.1999);
    CPPUNIT_ASSERT(!trigger.shouldWrite(t, tindex++));t += dt;
    CPPUNIT_ASSERT(trigger.shouldWrite(t, tindex++));t += dt;
    CPPUNIT_ASSERT(!trigger.shouldWrite(t, tindex++));t += dt;
    CPPUNIT_ASSERT(trigger.shouldWrite(t, tindex++));t += dt;
    CPPUNIT_ASSERT(!trigger.shouldWrite(t, tindex++));t += dt;

    trigger.setTimeSkip(0.2999);
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
