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

#include "PyreComponent.hh" // Implementation of class methods

#include "journals.hh"

#include <Python.h>

#include <stdexcept> // USES std::logic_error
#include <iostream> \
    // USES std::cerr

// ----------------------------------------------------------------------
// Constructor
pylith::utils::PyreComponent::PyreComponent(void) :
    _name(""),
    _identifier("unknown")
{ // constructor
    if (!Py_IsInitialized()) {
        throw std::logic_error("Python must be initialized to use PyreComponent in C++.");
    } // if
} // constructor


// ----------------------------------------------------------------------
// Destructor
pylith::utils::PyreComponent::~PyreComponent(void)
{ // destructor
} // destructor


// ----------------------------------------------------------------------
// Set name of component.
void
pylith::utils::PyreComponent::name(const char* value)
{ // name
    if (!strlen(value)) {
        throw std::logic_error("Cannot set name of Pyre component to empty string.");
    } // if
    _name = value;
} // name


// ----------------------------------------------------------------------
// Get name of component.
const char*
pylith::utils::PyreComponent::name(void) const
{ // name
    return _name.c_str();
} // name


// ----------------------------------------------------------------------
// Set component identifier (identifies object in component hierarchy).
void
pylith::utils::PyreComponent::identifier(const char* value)
{ // identifier
    if (!strlen(value)) {
        throw std::logic_error("Cannot set name of Pyre identifier to empty string.");
    } // if
    _identifier = value;
} // identifier


// ----------------------------------------------------------------------
// Get component identifier (identifies object in component hierarchy).
const char*
pylith::utils::PyreComponent::identifier(void) const
{ // identifier
    return _identifier.c_str();
} // identifier


// End of file
