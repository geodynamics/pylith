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

#include "ObservedSubject.hh" // Implementation of class methods

#include "pylith/feassemble/Observer.hh" // USES Observer
#include "pylith/topology/Field.hh" // USES Field

#include "pylith/utils/error.hh" // USES PYLITH_METHOD_BEGIN/END
#include "pylith/utils/journals.hh" // USES PYLITH_JOURNAL_DEBUG

#include <typeinfo> \
    // USES typeid()

// ----------------------------------------------------------------------
// Constructor.
pylith::feassemble::ObservedSubject::ObservedSubject(void) {}


// ----------------------------------------------------------------------
// Destructor
pylith::feassemble::ObservedSubject::~ObservedSubject(void) {
    deallocate();
} // destructor

// ----------------------------------------------------------------------
// Deallocate PETSc and local data structures.
void
pylith::feassemble::ObservedSubject::deallocate(void) {
    _observers.clear(); // Memory allocation of Observer* managed elsewhere.
} // deallocate


// ----------------------------------------------------------------------
// Register observer to receive notifications.
void
pylith::feassemble::ObservedSubject::registerObserver(pylith::feassemble::Observer* observer) {
    PYLITH_METHOD_BEGIN;
    PYLITH_JOURNAL_DEBUG("registerObserver(observer="<<typeid(observer).name()<<")");

    if (observer) {
        _observers.insert(observer);
    } // if

    PYLITH_METHOD_END;
} // registerObserver


// ----------------------------------------------------------------------
// Remove observer from receiving notifications.
void
pylith::feassemble::ObservedSubject::removeObserver(pylith::feassemble::Observer* observer) {
    PYLITH_METHOD_BEGIN;
    PYLITH_JOURNAL_DEBUG("removeObserver(observer="<<typeid(observer).name()<<")");

    if (observer) {
        _observers.erase(observer);
    } // if

    PYLITH_METHOD_END;
} // removeObserver


// ----------------------------------------------------------------------
// Notify observers.
void
pylith::feassemble::ObservedSubject::notifyObservers(const PylithReal t,
                                                     const PylithInt tindex,
                                                     const pylith::topology::Field& solution,
                                                     const bool infoOnly) {
    PYLITH_METHOD_BEGIN;
    PYLITH_JOURNAL_DEBUG("notifyObservers(t="<<t<<", tindex="<<tindex<<", solution="<<solution.label()<<")");

    typedef std::set<Observer*>::iterator observers_iter;
    for (observers_iter iter = _observers.begin(); iter != _observers.end(); ++iter) {
        assert(*iter);
        (*iter)->update(t, tindex, solution);
    } // for

    PYLITH_METHOD_END;
} // notifyObservers



// End of file
