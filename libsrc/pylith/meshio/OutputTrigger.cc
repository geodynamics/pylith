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
#include <stdexcept>

#include "OutputTrigger.hh" // Implementation of class methods

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
        std::ostringstream msg;
        msg << "Time scale ("<<value<<") for solution observer is nonpositive.";
        throw std::logic_error(msg.str());
    } // if
    _timeScale = value;
} // setTimeScale


// End of file
