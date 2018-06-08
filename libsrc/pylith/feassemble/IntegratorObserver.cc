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

#include "IntegratorObserver.hh" // Implementation of class methods

#include "pylith/feassemble/IntegratorPointwise.hh" \
    // USES IntegratorPointwise


// ----------------------------------------------------------------------
// Constructor with integrator to observe.
pylith::feassemble::IntegratorObserver::IntegratorObserver(pylith::feassemble::IntegratorPointwise* const integrator) :
    _integrator(integrator)
{ // constructor
    if (_integrator) {
        _integrator->registerObserver(this);
    } // if
} // constructor


// ----------------------------------------------------------------------
// Destructor
pylith::feassemble::IntegratorObserver::~IntegratorObserver(void) {}


// ----------------------------------------------------------------------
// Deallocate PETSc and local data structures.
void
pylith::feassemble::IntegratorObserver::deallocate(void) {
    Observer::deallocate();

    if (_integrator) {
        _integrator->removeObserver(this);
    } // if
} // deallocate


// End of file
