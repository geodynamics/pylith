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

#include "ObserversPhysics.hh" // Implementation of class methods

#include "pylith/problems/ObserverPhysics.hh" // USES ObserverPhysics
#include "pylith/feassemble/PhysicsImplementation.hh" // USES PhysicsImplementation
#include "pylith/topology/Field.hh" // USES Field

#include "pylith/utils/error.hh" // USES PYLITH_METHOD_BEGIN/END
#include "pylith/utils/journals.hh" // USES PYLITH_JOURNAL_DEBUG

#include <typeinfo> // USES typeid()

// ------------------------------------------------------------------------------------------------
// Constructor.
pylith::problems::ObserversPhysics::ObserversPhysics(void) {
    // GenericComponent::setName("observersphysics");
} // constructor


// ------------------------------------------------------------------------------------------------
// Destructor
pylith::problems::ObserversPhysics::~ObserversPhysics(void) {
    deallocate();
} // destructor


// ------------------------------------------------------------------------------------------------
// Deallocate PETSc and local data structures.
void
pylith::problems::ObserversPhysics::deallocate(void) {
    _observers.clear(); // Memory allocation of Observer* managed elsewhere.
} // deallocate


// ------------------------------------------------------------------------------------------------
// Register observer to receive notifications.
void
pylith::problems::ObserversPhysics::registerObserver(pylith::problems::ObserverPhysics* observer) {
    PYLITH_METHOD_BEGIN;
    PYLITH_JOURNAL_DEBUG("registerObserver(observer="<<typeid(observer).name()<<")");

    if (observer) {
        _observers.insert(observer);
    } // if

    PYLITH_METHOD_END;
} // registerObserver


// ------------------------------------------------------------------------------------------------
// Remove observer from receiving notifications.
void
pylith::problems::ObserversPhysics::removeObserver(pylith::problems::ObserverPhysics* observer) {
    PYLITH_METHOD_BEGIN;
    PYLITH_JOURNAL_DEBUG("removeObserver(observer="<<typeid(observer).name()<<")");

    if (observer) {
        _observers.erase(observer);
    } // if

    PYLITH_METHOD_END;
} // removeObserver


// ------------------------------------------------------------------------------------------------
// Get number of observers.
size_t
pylith::problems::ObserversPhysics::size(void) const {
    return _observers.size();
} // count`


// ------------------------------------------------------------------------------------------------
// Set physics implementation in observers (for callbacks)
void
pylith::problems::ObserversPhysics::setPhysicsImplementation(const pylith::feassemble::PhysicsImplementation* const physics) {
    PYLITH_METHOD_BEGIN;
    PYLITH_JOURNAL_DEBUG("setPhysicsImplementation(physics="<<physics<<")");

    for (iterator iter = _observers.begin(); iter != _observers.end(); ++iter) {
        assert(*iter);
        (*iter)->setPhysicsImplementation(physics);
    } // for

    PYLITH_METHOD_END;
} // setPhysicsImplemetation


// ------------------------------------------------------------------------------------------------
// Set time scale in observers.
void
pylith::problems::ObserversPhysics::setTimeScale(const PylithReal value) {
    PYLITH_METHOD_BEGIN;
    PYLITH_JOURNAL_DEBUG("setTimeScale(value="<<value<<")");

    for (iterator iter = _observers.begin(); iter != _observers.end(); ++iter) {
        assert(*iter);
        (*iter)->setTimeScale(value);
    } // for

    PYLITH_METHOD_END;
} // setTimeScale


// ------------------------------------------------------------------------------------------------
// Verify observers.
void
pylith::problems::ObserversPhysics::verifyObservers(const pylith::topology::Field& solution) const {
    PYLITH_METHOD_BEGIN;
    PYLITH_JOURNAL_DEBUG("verifyObservers(solution="<<solution.getLabel()<<")");

    for (iterator iter = _observers.begin(); iter != _observers.end(); ++iter) {
        assert(*iter);
        (*iter)->verifyConfiguration(solution);
    } // for

    PYLITH_METHOD_END;
} // verifyObservers


// ------------------------------------------------------------------------------------------------
// Notify observers.
void
pylith::problems::ObserversPhysics::notifyObservers(const PylithReal t,
                                                    const PylithInt tindex,
                                                    const pylith::topology::Field& solution,
                                                    const bool infoOnly) {
    PYLITH_METHOD_BEGIN;
    PYLITH_JOURNAL_DEBUG("notifyObservers(t="<<t<<", tindex="<<tindex<<", solution="<<solution.getLabel()<<")");

    for (iterator iter = _observers.begin(); iter != _observers.end(); ++iter) {
        assert(*iter);
        (*iter)->update(t, tindex, solution, infoOnly);
    } // for

    PYLITH_METHOD_END;
} // notifyObservers


// End of file
