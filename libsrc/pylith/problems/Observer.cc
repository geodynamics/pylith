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

#include "Observer.hh" \
    // Implementation of class methods


// ----------------------------------------------------------------------
// Constructor.
pylith::problems::Observer::Observer(void) {}


// ----------------------------------------------------------------------
// Destructor
pylith::problems::Observer::~Observer(void) {
    deallocate();
} // destructor


// ----------------------------------------------------------------------
// Deallocate PETSc and local data structures.
void
pylith::problems::Observer::deallocate(void) {}


// End of file
