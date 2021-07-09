// -*- C++ -*-
//
// ======================================================================
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
// ======================================================================
//

#include <portinfo>

#include "OutputTriggerTime.hh" // Implementation of class methods

#include "pylith/utils/constdefs.h" // USES PYLITH_MAXSCALAR
#include "pylith/utils/error.hh" // USES PYLITH_METHOD_BEGIN/END
#include "pylith/utils/journals.hh" // USES PYLITH_COMPONENT_*

// ---------------------------------------------------------------------------------------------------------------------
// Constructor
pylith::meshio::OutputTriggerTime::OutputTriggerTime(void) :
    _timeSkip(0.0),
    _timeNondimWrote(-PYLITH_MAXSCALAR) {
    PyreComponent::setName("outputtriggertime");
} // constructor


// ---------------------------------------------------------------------------------------------------------------------
// Destructor
pylith::meshio::OutputTriggerTime::~OutputTriggerTime(void) {}


// ---------------------------------------------------------------------------------------------------------------------
// Set elapsed time between writes.
void
pylith::meshio::OutputTriggerTime::setTimeSkip(const double value) {
    PYLITH_COMPONENT_DEBUG("OutputTriggerTime::setTimeSkip(value="<<value<<")");

    _timeSkip = (value >= 0.0) ? value : 0.0;
} // setTimeSkip


// ---------------------------------------------------------------------------------------------------------------------
// Get elapsed time between writes.
double
pylith::meshio::OutputTriggerTime::getTimeSkip(void) const {
    return _timeSkip;
} // getTimeSkip


// ---------------------------------------------------------------------------------------------------------------------
// Check whether we want to write output at time t.
bool
pylith::meshio::OutputTriggerTime::shouldWrite(const PylithReal t,
                                               const PylithInt timeStep) {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("OutputTriggerTime::shouldWrite(t="<<t<<", timeStep="<<timeStep<<")");

    bool isWrite = false;
    if (t - _timeNondimWrote >= _timeSkip / _timeScale) {
        isWrite = true;
        _timeNondimWrote = t;
    } // if

    PYLITH_METHOD_RETURN(isWrite);
} // shouldWrite


// End of file
