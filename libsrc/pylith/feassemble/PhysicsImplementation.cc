// -*- C++ -*-
//
// ----------------------------------------------------------------------
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
// ----------------------------------------------------------------------
//

#include <portinfo>

#include "pylith/feassemble/PhysicsImplementation.hh" // implementation of object methods

#include "pylith/topology/Mesh.hh" // USES Mesh
#include "pylith/topology/Field.hh" // USES Field
#include "pylith/problems/ObserversPhysics.hh" // USES ObserversPhysics
#include "pylith/problems/Physics.hh" // USES Physics

#include "pylith/utils/EventLogger.hh" // USES EventLogger
#include "pylith/utils/error.hh" // USES PYLITH_METHOD_*
#include "pylith/utils/journals.hh" // USES PYLITH_JOURNAL_*

#include <cassert> // USES assert()
#include <stdexcept> // USES std::runtime_error

// ------------------------------------------------------------------------------------------------
// Default constructor.
pylith::feassemble::PhysicsImplementation::PhysicsImplementation(pylith::problems::Physics* const physics) :
    _physics(physics),
    _auxiliaryField(NULL),
    _derivedField(NULL),
    _observers(NULL),
    _logger(NULL)
{}


// ------------------------------------------------------------------------------------------------
// Destructor.
pylith::feassemble::PhysicsImplementation::~PhysicsImplementation(void) {
    deallocate();
} // destructor


// ------------------------------------------------------------------------------------------------
// Deallocate PETSc and local data structures.
void
pylith::feassemble::PhysicsImplementation::deallocate(void) {
    PYLITH_METHOD_BEGIN;

    delete _auxiliaryField;_auxiliaryField = NULL;
    delete _derivedField;_derivedField = NULL;
    _observers = NULL; // :KLUDGE: Use shared pointer. _observers held by physics.
    delete _logger;_logger = NULL;

    PYLITH_METHOD_END;
} // deallocate


// ------------------------------------------------------------------------------------------------
// Get auxiliary field.
const pylith::topology::Field*
pylith::feassemble::PhysicsImplementation::getAuxiliaryField(void) const {
    return _auxiliaryField;
} // getAuxiliaryField


// ------------------------------------------------------------------------------------------------
// Get derived field.
const pylith::topology::Field*
pylith::feassemble::PhysicsImplementation::getDerivedField(void) const {
    return _derivedField;
} // getDerivedField


// ------------------------------------------------------------------------------------------------
// Notify observers of current solution.
void
pylith::feassemble::PhysicsImplementation::notifyObservers(const PylithReal t,
                                                           const PylithInt tindex,
                                                           const pylith::topology::Field& solution) {
    if (!_observers) {
        return;
    } // if

    assert(_observers);
    const bool infoOnly = false;
    _observers->notifyObservers(t, tindex, solution, infoOnly);
} // _notifyObservers


// End of file
