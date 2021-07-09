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

#include "PyreComponent.hh" // Implementation of class methods

#include "journals.hh"

#include <Python.h>

#include <stdexcept> // USES std::logic_error
#include <iostream> // USES std::cerr

// ----------------------------------------------------------------------
// Constructor
pylith::utils::PyreComponent::PyreComponent(void) :
    _name(""),
    _identifier("unknown") { // constructor
    if (!Py_IsInitialized()) {
        throw std::logic_error("Python must be initialized to use PyreComponent in C++.");
    } // if
} // constructor


// ----------------------------------------------------------------------
// Destructor
pylith::utils::PyreComponent::~PyreComponent(void) {}


// ----------------------------------------------------------------------
// Set name of component.
void
pylith::utils::PyreComponent::setName(const char* value) {
    if (!strlen(value)) {
        throw std::logic_error("Cannot set name of Pyre component to empty string.");
    } // if
    _name = value;
} // setName


// ----------------------------------------------------------------------
// Get name of component.
const char*
pylith::utils::PyreComponent::getName(void) const {
    return _name.c_str();
} // getName


// ----------------------------------------------------------------------
// Set component identifier (identifies object in component hierarchy).
void
pylith::utils::PyreComponent::setIdentifier(const char* value) {
    if (!strlen(value)) {
        throw std::logic_error("Cannot set name of Pyre identifier to empty string.");
    } // if
    _identifier = value;
} // setIdentifier


// ----------------------------------------------------------------------
// Get component identifier (identifies object in component hierarchy).
const char*
pylith::utils::PyreComponent::getIdentifier(void) const {
    return _identifier.c_str();
} // getIdentifier


// End of file
