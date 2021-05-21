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
// Copyright (c) 2010-2017 University of California, Davis
//
// See COPYING for license information.
//
// ======================================================================
//

#include <portinfo>

#include "Fields.hh" // Implementation of class methods

#include "Field.hh" // USES Field

#include <pylith/utils/error.hh> // USES PYLITH_CHECK_ERROR

#include <sstream> // USES std::ostringstream
#include <stdexcept> // USES std::runtime_error

// ----------------------------------------------------------------------
// Default constructor.
pylith::topology::Fields::Fields(const Mesh& mesh) :
    _mesh(mesh) { // constructor
} // constructor


// ----------------------------------------------------------------------
// Destructor.
pylith::topology::Fields::~Fields(void) {
    deallocate();
} // destructor


// ----------------------------------------------------------------------
// Deallocate PETSc and local data structures.
void
pylith::topology::Fields::deallocate(void) {
    PYLITH_METHOD_BEGIN;

    const map_type::iterator begin = _fields.begin();
    const map_type::iterator end = _fields.end();
    for (map_type::iterator iter = begin; iter != end; ++iter) {
        delete iter->second;iter->second = NULL;
    } // for
    _fields.clear();

    PYLITH_METHOD_END;
} // deallocate


// ----------------------------------------------------------------------
// Check if fields contains a given field.
bool
pylith::topology::Fields::hasField(const char* name) const {
    PYLITH_METHOD_BEGIN;
    assert(name);

    map_type::const_iterator iter = _fields.find(name);

    PYLITH_METHOD_RETURN(iter != _fields.end());
} // hasField


// ----------------------------------------------------------------------
// Add field.
void
pylith::topology::Fields::add(const char* name,
                              const char* label) {
    PYLITH_METHOD_BEGIN;
    assert(name);
    assert(label);

    if (hasField(name)) {
        std::ostringstream msg;
        msg << "Could not add field '" << name << "' to fields manager, because it already exists.";
        throw std::runtime_error(msg.str());
    } // if

    _fields[name] = new Field(_mesh);
    _fields[name]->setLabel(label);

    PYLITH_METHOD_END;
} // add


// ----------------------------------------------------------------------
// Delete field.
void
pylith::topology::Fields::del(const char* name) {
    PYLITH_METHOD_BEGIN;
    assert(name);

    map_type::iterator iter = _fields.find(name);
    if (iter == _fields.end()) {
        std::ostringstream msg;
        msg << "Could not find field '" << name << "' in fields manager to delete.";
        throw std::runtime_error(msg.str());
    } // if
    delete iter->second;iter->second = 0;
    _fields.erase(name);

    PYLITH_METHOD_END;
} // del


// ----------------------------------------------------------------------
// Get field.
const pylith::topology::Field&
pylith::topology::Fields::get(const char* name) const {
    PYLITH_METHOD_BEGIN;
    assert(name);

    map_type::const_iterator iter = _fields.find(name);
    if (iter == _fields.end()) {
        std::ostringstream msg;
        msg << "Could not find field '" << name << "' in fields manager for retrieval.";
        throw std::runtime_error(msg.str());
    } // if

    PYLITH_METHOD_RETURN(*iter->second);
} // get


// ----------------------------------------------------------------------
// Get field.
pylith::topology::Field&
pylith::topology::Fields::get(const char* name) {
    PYLITH_METHOD_BEGIN;
    assert(name);

    map_type::iterator iter = _fields.find(name);
    if (iter == _fields.end()) {
        std::ostringstream msg;
        msg << "Could not find field '" << name << "' in fields manager for retrieval.";
        throw std::runtime_error(msg.str());
    } // if

    PYLITH_METHOD_RETURN(*iter->second);
} // get


// ----------------------------------------------------------------------
// Get mesh associated with fields.
const pylith::topology::Mesh&
pylith::topology::Fields::mesh(void) const {
    return _mesh;
} // mesh


// End of file
