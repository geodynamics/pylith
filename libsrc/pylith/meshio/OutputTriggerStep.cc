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

#include "OutputTriggerStep.hh" // Implementation of class methods

#include "pylith/utils/constdefs.h" // USES PYLITH_MAXSCALAR
#include "pylith/utils/error.hh" // USES PYLITH_METHOD_BEGIN/END
#include "pylith/utils/journals.hh" // USES PYLITH_COMPONENT_*

// ---------------------------------------------------------------------------------------------------------------------
// Constructor
pylith::meshio::OutputTriggerStep::OutputTriggerStep(void) :
    _numStepsSkip(0),
    _stepWrote(PYLITH_MININT+10) { // constructor
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
