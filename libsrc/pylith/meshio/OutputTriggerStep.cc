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

#include "OutputTriggerStep.hh" // Implementation of class methods

#include "pylith/utils/constdefs.h" // USES PYLITH_MAXSCALAR
#include "pylith/utils/error.hh" // USES PYLITH_METHOD_BEGIN/END
#include "pylith/utils/journals.hh" // USES PYLITH_COMPONENT_*

// ----------------------------------------------------------------------
const char* pylith::meshio::OutputTriggerStep::_pyreComponent = "outputtriggerstep";

// ----------------------------------------------------------------------
// Constructor
pylith::meshio::OutputTriggerStep::OutputTriggerStep(void) :
    _numTimeStepsSkip(0),
    _timeStepWrote(PYLITH_MININT+10)
{ // constructor
    PyreComponent::name(_pyreComponent);
} // constructor

// ----------------------------------------------------------------------
// Destructor
pylith::meshio::OutputTriggerStep::~OutputTriggerStep(void) {}

// ----------------------------------------------------------------------
// Set number of time steps to skip between writes.
void
pylith::meshio::OutputTriggerStep::numTimeStepsSkip(const int value) {
    PYLITH_COMPONENT_DEBUG("OutputTriggerStep::numTimeStepsSkip(value="<<value<<")");

    _numTimeStepsSkip = (value >= 0) ? value : 0;
} // numTimeStepsSkip

// ----------------------------------------------------------------------
// Get number of time steps to skip between writes.
int
pylith::meshio::OutputTriggerStep::numTimeStepsSkip(void) const {
    return _numTimeStepsSkip;
} // numTimeStepsSkip

// ----------------------------------------------------------------------
// Check whether we want to write output at time t.
bool
pylith::meshio::OutputTriggerStep::shouldWrite(const PylithReal t,
                                               const PylithInt timeStep) {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("OutputTriggerStep::shouldWrite(t="<<t<<", timeStep="<<timeStep<<")");

    bool shouldWrite = false;
    if (timeStep - _timeStepWrote > _numTimeStepsSkip) {
        shouldWrite = true;
        _timeStepWrote = timeStep;
    } // if

    PYLITH_METHOD_RETURN(shouldWrite);
} // shouldWrite


// End of file
