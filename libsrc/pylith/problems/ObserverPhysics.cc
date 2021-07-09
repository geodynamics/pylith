// -*- C++ -*-
//
// ======================================================================
//
// Brad T. Aagaard, U.S. Geological Survey
// Charles A. Williams, GNS Science
// Matthew G. Knepley, University at Buffalo
//
// This code was developed as part of the Computational Infrastructure
// for Geodynamics (http://geodynamics.org).
//
// Copyright (c) 2010-2021 University of California, Davis
//
// See LICENSE.md for license information.
//
// ======================================================================
//

#include <portinfo>

#include "ObserverPhysics.hh" // Implementation of class methods

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
