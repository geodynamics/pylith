// =================================================================================================
// This code is part of PyLith, developed through the Computational Infrastructure
// for Geodynamics (https://github.com/geodynamics/pylith).
//
// Copyright (c) 2010-2024, University of California, Davis and the PyLith Development Team.
// All rights reserved.
//
// See https://mit-license.org/ and LICENSE.md and for license information.
// =================================================================================================

#include <portinfo>

#include "pylith/utils/GenericComponent.hh" // Implementation of class methods

#include "pylith/utils/journals.hh"

#include <Python.h>

#include <stdexcept> // USES std::logic_error
#include <iostream>

// ----------------------------------------------------------------------
// Constructor
pylith::utils::GenericComponent::GenericComponent(void) :
    _name("") {
    if (!Py_IsInitialized()) {
        throw std::logic_error("Python must be initialized to use GenericComponent in C++.");
    } // if
} // constructor


// ----------------------------------------------------------------------
// Destructor
pylith::utils::GenericComponent::~GenericComponent(void) {}


// ----------------------------------------------------------------------
// Set name of component.
void
pylith::utils::GenericComponent::setName(const char* value) {
    if (!strlen(value)) {
        throw std::logic_error("Cannot set name of Generic component to empty string.");
    } // if
    _name = value;
} // setName


// ----------------------------------------------------------------------
// Get name of component.
const char*
pylith::utils::GenericComponent::getName(void) const {
    return _name.c_str();
} // getName


// End of file
