// =================================================================================================
// This code is part of PyLith, developed through the Computational Infrastructure
// for Geodynamics (https://github.com/geodynamics/pylith).
//
// Copyright (c) 2010-2025, University of California, Davis and the PyLith Development Team.
// All rights reserved.
//
// See https://mit-license.org/ and LICENSE.md and for license information.
// =================================================================================================

#include <portinfo>

#include "pylith/initializers/InitializePhase.hh" // implementation of class methods

// ------------------------------------------------------------------------------------------------
// Default constructor
pylith::initializers::InitializePhase::InitializePhase(void) {}


// ------------------------------------------------------------------------------------------------
// Default destructor
pylith::initializers::InitializePhase::~InitializePhase(void) {
    deallocate();
} // destructor


// ------------------------------------------------------------------------------------------------
// Deallocate PETSc and local data structures.
void
pylith::initializers::InitializePhase::deallocate(void) {
}


// End of file
