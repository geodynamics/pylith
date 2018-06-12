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

#include "PetscComponent.hh" \
    // Implementation of class methods

// ----------------------------------------------------------------------
// Constructor
pylith::utils::PetscComponent::PetscComponent(void) {}

// ----------------------------------------------------------------------
// Destructor
pylith::utils::PetscComponent::~PetscComponent(void) {
    deallocate();
} // destructor


// End of file
