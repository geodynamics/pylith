// =================================================================================================
// This code is part of PyLith, developed through the Computational Infrastructure
// for Geodynamics (https://github.com/geodynamics/pylith).
//
// Copyright (c) 2010-2026, University of California, Davis and the PyLith Development Team.
// All rights reserved.
//
// See https://mit-license.org/ and LICENSE.md and for license information.
// =================================================================================================
#include <portinfo>

#include "pylith/meshio/OutputTrigger.hh" // Implementation of class methods

#include "pylith/utils/journals.hh" // USES journal macros
#include "pylith/utils/Exceptions.hh" // USES Exception


// ---------------------------------------------------------------------------------------------------------------------
// Constructor
pylith::meshio::OutputTrigger::OutputTrigger(void) :
    _timeScale(1.0) {}


// ---------------------------------------------------------------------------------------------------------------------
// Destructor
pylith::meshio::OutputTrigger::~OutputTrigger(void) {}


// ---------------------------------------------------------------------------------------------------------------------
// Set time scale.
void
pylith::meshio::OutputTrigger::setTimeScale(const PylithReal value) {
    if (value <= 0.0) {
        PYLITH_COMPONENT_ERROR(pylith::ValueError, pylith::journal::user_input,
                               "Time scale ("<<value<<") for solution observer is nonpositive.");
    } // if
    _timeScale = value;
} // setTimeScale


// End of file
