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
// Copyright (c) 2010-2015 University of California, Davis
//
// See COPYING for license information.
//
// ======================================================================
//

#include <portinfo>

#include "InitialConditions.hh" // implementation of class methods

#include "pylith/utils/error.hh" // USES PYLITH_CHECK_ERROR

// ----------------------------------------------------------------------
// Constructor
pylith::problems::InitialConditions::InitialConditions(void)
{}


// ---------------------------------------------------------------------------------------------------------------------
// Destructor
pylith::problems::InitialConditions::~InitialConditions(void) {
    deallocate();
} // destructor


// ---------------------------------------------------------------------------------------------------------------------
// Deallocate PETSc and local data structures.
void
pylith::problems::InitialConditions::deallocate(void) {
    PYLITH_METHOD_BEGIN;

    PYLITH_METHOD_END;
} // deallocate


// End of file
