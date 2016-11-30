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

#include "journal/debug.h"
#include "journal/info.h"
#include "journal/warning.h"
#include "journal/error.h"

#include <cassert> // USES assert()
#include <stdexcept> // USES std::logic_error

// ----------------------------------------------------------------------
// Constructor
pylith::utils::JournalingComponent::JournalingComponent(void) :
    _debug(0),
    _info(0),
    _warning(0),
    _error(0),
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
    delete _debug; _debug = 0;
    delete _info; _info = 0;
    delete _warning; _warning = 0;
    delete _error; _error = 0;
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
/** Get name of journal.
 *
 * @returns Name of journal.
 */
const char*
pylith::utils::JournalingComponent::name(void) const
{ // name
    return _name.c_str();
} // name

// ----------------------------------------------------------------------
// Setup journals.
void
pylith::utils::JournalingComponent::initialize(void)
{ // initialize
    if (!_name.length()) {
        throw std::logic_error("Name of journaling component not set.");
    } // if

    delete _debug; _debug = new journal::debug_t(_name); assert(_debug);
    delete _info; _info = new journal::info_t(_name); assert(_info);
    delete _warning; _warning = new journal::warning_t(_name); assert(_warning);
    delete _error; _error = new journal::error_t(_name); assert(_error);
} // initialize

// ----------------------------------------------------------------------
// Get debug journal.
journal::debug_t&
pylith::utils::JournalingComponent::debug(void)
{ // debug
    assert(_debug);
    return *_debug;
} // debug

// ----------------------------------------------------------------------
// Get info journal.
journal::info_t&
pylith::utils::JournalingComponent::info(void)
{ // info
    assert(_info);
    return *_info;
} // info

// ----------------------------------------------------------------------
// Get warning journal.
journal::warning_t&
pylith::utils::JournalingComponent::warning(void)
{ // warning
    assert(_warning);
    return *_warning;
} // warning

// ----------------------------------------------------------------------
// Get error journal.
journal::error_t&
pylith::utils::JournalingComponent::error(void)
{ // error
    assert(_error);
    return *_error;
} // error


// End of file
