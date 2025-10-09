// =================================================================================================
// This code is part of PyLith, developed through the Computational Infrastructure
// for Geodynamics (https://github.com/geodynamics/pylith).
//
// Copyright (c) 2010-2025, University of California, Davis and the PyLith Development Team.
// All rights reserved.
//
// See https://mit-license.org/ and LICENSE.md and for license information.
// =================================================================================================

#include <portinfo>

#include "pylith/meshio/OutputTriggerStep.hh" // Implementation of class methods

#include "pylith/utils/constants.hh" // USES pylith::min_int
#include "pylith/utils/error.hh" // USES PYLITH_METHOD_BEGIN/END
#include "pylith/utils/journals.hh" // USES PYLITH_COMPONENT_*

// ---------------------------------------------------------------------------------------------------------------------
// Constructor
pylith::meshio::OutputTriggerStep::OutputTriggerStep(void) :
    _numStepsSkip(0),
    _stepWrote(pylith::min_int+10) { // constructor
    PyreComponent::setName("outputtriggerstep");
} // constructor


// ---------------------------------------------------------------------------------------------------------------------
// Destructor
pylith::meshio::OutputTriggerStep::~OutputTriggerStep(void) {}


// ---------------------------------------------------------------------------------------------------------------------
// Set number of steps to skip between writes.
void
pylith::meshio::OutputTriggerStep::setNumStepsSkip(const int value) {
    PYLITH_COMPONENT_DEBUG("OutputTriggerStep::setNumStepsSkip(value="<<value<<")");

    _numStepsSkip = (value >= 0) ? value : 0;
} // setNumStepsSkip


// ---------------------------------------------------------------------------------------------------------------------
// Get number of steps to skip between writes.
int
pylith::meshio::OutputTriggerStep::getNumStepsSkip(void) const {
    return _numStepsSkip;
} // getNumStepsSkip


// ---------------------------------------------------------------------------------------------------------------------
// Check whether we want to write output at time t.
bool
pylith::meshio::OutputTriggerStep::shouldWrite(const PylithReal t,
                                               const PylithInt tindex) {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("OutputTriggerStep::shouldWrite(t="<<t<<", tindex="<<tindex<<")");

    bool isWrite = false;
    if (tindex - _stepWrote > _numStepsSkip) {
        isWrite = true;
        _stepWrote = tindex;
    } // if

    PYLITH_METHOD_RETURN(isWrite);
} // shouldWrite


// End of file
