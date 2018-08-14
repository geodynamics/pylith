// -*- C++ -*-
//
// ======================================================================
//
// Brad T. Aagaard, U.S. Geological Survey
// Charles A. Williams, GNS Science
// Matthew G. Knepley, University of Chicago
//
// This code was developed as part of the Computational Infrastructure
// for Geodynamics (http://geodynamics.org).
//
// Copyright (c) 2010-2016 University of California, Davis
//
// See COPYING for license information.
//
// ======================================================================
//

#include <portinfo>

#include "OutputTriggerTime.hh" // Implementation of class methods

#include "pylith/utils/constdefs.h" // USES PYLITH_MAXSCALAR
#include "pylith/utils/error.hh" // USES PYLITH_METHOD_BEGIN/END
#include "pylith/utils/journals.hh" // USES PYLITH_COMPONENT_*

// ----------------------------------------------------------------------
const char* pylith::meshio::OutputTriggerTime::_pyreComponent = "outputtriggertime";

// ----------------------------------------------------------------------
// Constructor
pylith::meshio::OutputTriggerTime::OutputTriggerTime(void) :
    _timeSkip(0.0),
    _timeWrote(-PYLITH_MAXSCALAR)
{ // constructor
    PyreComponent::name(_pyreComponent);
} // constructor

// ----------------------------------------------------------------------
// Destructor
pylith::meshio::OutputTriggerTime::~OutputTriggerTime(void) {}

// ----------------------------------------------------------------------
// Set elapsed time between writes.
void
pylith::meshio::OutputTriggerTime::timeSkip(const double value) {
    PYLITH_COMPONENT_DEBUG("OutputTriggerTime::timeSkip(value="<<value<<")");

    _timeSkip = (value >= 0.0) ? value : 0.0;
} // timeSkip

// ----------------------------------------------------------------------
// Get elapsed time between writes.
double
pylith::meshio::OutputTriggerTime::timeSkip(void) const {
    return _timeSkip;
} // timeSkip


// ----------------------------------------------------------------------
// Check whether we want to write output at time t.
bool
pylith::meshio::OutputTriggerTime::shouldWrite(const PylithReal t,
                                               const PylithInt timeStep) {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("OutputTriggerTime::shouldWrite(t="<<t<<", timeStep="<<timeStep<<")");

    bool shouldWrite = false;
    if (t - _timeWrote >= _timeSkip) {
        shouldWrite = true;
        _timeWrote = t;
    } // if

    PYLITH_METHOD_RETURN(shouldWrite);
} // shouldWrite


// End of file
