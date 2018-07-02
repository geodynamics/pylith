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

#include "FieldFilterNone.hh" \
    // Implementation of class methods

// ----------------------------------------------------------------------
// Constructor
pylith::meshio::FieldFilterNone::FieldFilterNone(void) {}

// ----------------------------------------------------------------------
// Destructor
pylith::meshio::FieldFilterNone::~FieldFilterNone(void) {}

// ----------------------------------------------------------------------
// Copy constructor.
pylith::meshio::FieldFilterNone::FieldFilterNone(const FieldFilterNone& f) :
    FieldFilter(f)
{}

// ----------------------------------------------------------------------
// Create copy of filter.
pylith::meshio::FieldFilter*
pylith::meshio::FieldFilterNone::clone(void) const {
    return new FieldFilterNone(*this);
} // clone

// ----------------------------------------------------------------------
// Filter field.
pylith::topology::Field*
pylith::meshio::FieldFilterNone::filter(pylith::topology::Field* fieldIn) {
    return fieldIn;
} // filter


// End of file
