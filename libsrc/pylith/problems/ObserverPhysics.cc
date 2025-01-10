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

#include "pylith/problems/ObserverPhysics.hh" // Implementation of class methods

#include "pylith/utils/error.hh" // USES PYLITH_METHOD_BEGIN/END

// ------------------------------------------------------------------------------------------------
// Constructor.
pylith::problems::ObserverPhysics::ObserverPhysics(void) :
    _physics(NULL) {}


// ------------------------------------------------------------------------------------------------
// Destructor
pylith::problems::ObserverPhysics::~ObserverPhysics(void) {
    deallocate();
} // destructor


// ------------------------------------------------------------------------------------------------
// Deallocate PETSc and local data structures.
void
pylith::problems::ObserverPhysics::deallocate(void) {
    _physics = NULL; // :TODO: Use shared pointer.
} // deallocate


// ------------------------------------------------------------------------------------------------
// Set physics implementation to observe.
void
pylith::problems::ObserverPhysics::setPhysicsImplementation(const pylith::feassemble::PhysicsImplementation* const physics) {
    PYLITH_METHOD_BEGIN;

    _physics = physics;

    PYLITH_METHOD_END;
} // setPhysicsImplemetation


// End of file
