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

#include "JournalingComponent.hh" // Implementation of class methods

#include <Python.h>

#include <stdexcept> // USES std::logic_error

// ----------------------------------------------------------------------
// Constructor
pylith::utils::JournalingComponent::JournalingComponent(void) :
    _name("")
{ // constructor
    if (!Py_IsInitialized()) {
        throw std::logic_error("Python must be initialized to use Pyre journals in C++.");
    } // if
} // constructor

// ----------------------------------------------------------------------
// Destructor
pylith::utils::JournalingComponent::~JournalingComponent(void)
{ // destructor
} // destructor

// ----------------------------------------------------------------------
// Set name of journal.
void
pylith::utils::JournalingComponent::name(const char* value)
{ // name
    if (!strlen(value)) {
        throw std::logic_error("Cannot set name of journaling component to empty string.");
    } // if
    _name = value;
} // name

// ----------------------------------------------------------------------
// Get name of journal.
const char*
pylith::utils::JournalingComponent::name(void) const
{ // name
    return _name.c_str();
} // name


// End of file
