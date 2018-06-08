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

#include "ProblemObserver.hh" // Implementation of class methods

#include "pylith/problems/Problem.hh" \
    // USES Problem


// ----------------------------------------------------------------------
// Constructor with problem to observe.
pylith::problems::ProblemObserver::ProblemObserver(pylith::problems::Problem* const problem) :
    _problem(problem)
{ // constructor
    if (_problem) {
        _problem->registerObserver(this);
    } // if
} // constructor


// ----------------------------------------------------------------------
// Destructor
pylith::problems::ProblemObserver::~ProblemObserver(void) {}


// ----------------------------------------------------------------------
// Deallocate PETSc and local data structures.
void
pylith::problems::ProblemObserver::deallocate(void) {
    Observer::deallocate();

    if (_problem) {
        _problem->removeObserver(this);
    } // if
} // deallocate


// End of file
