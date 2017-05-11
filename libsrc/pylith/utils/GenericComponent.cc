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

#include "GenericComponent.hh" // Implementation of class methods

#include "journals.hh"

#include <Python.h>

#include <stdexcept> // USES std::logic_error
#include <iostream>

// ----------------------------------------------------------------------
// Constructor
pylith::utils::GenericComponent::GenericComponent(void) :
    _name("")
{ // constructor
    if (!Py_IsInitialized()) {
        throw std::logic_error("Python must be initialized to use GenericComponent in C++.");
    } // if
} // constructor

// ----------------------------------------------------------------------
// Destructor
pylith::utils::GenericComponent::~GenericComponent(void)
{ // destructor
} // destructor


// ----------------------------------------------------------------------
// Set name of component.
void
pylith::utils::GenericComponent::name(const char* value)
{ // name
    if (!strlen(value)) {
        throw std::logic_error("Cannot set name of Generic component to empty string.");
    } // if
    _name = value;
} // name


// ----------------------------------------------------------------------
// Get name of component.
const char*
pylith::utils::GenericComponent::name(void) const
{ // name
    return _name.c_str();
} // name


// End of file
