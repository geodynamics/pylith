// -*- C++ -*-
//
// ----------------------------------------------------------------------
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
// ----------------------------------------------------------------------
//

#include <portinfo>

#include "DirichletTimeDependent.hh" // implementation of object methods

#include "pylith/topology/Field.hh" // USES Field

#include "pylith/utils/error.hh" // USES PYLITH_METHOD_BEGIN/END
#include "pylith/utils/journals.hh" // USES PYLITH_JOURNAL_*

#include <cassert> // USES assert()
#include <stdexcept> // USES std::runtime_error
#include <sstream> // USES std::ostringstream

// ----------------------------------------------------------------------
// Default constructor.
pylith::bc::DirichletTimeDependent::DirichletTimeDependent(void) :
    _useInitial(true),
    _useRate(false),
    _useTimeHistory(false)
{ // constructor
    JournalingComponent::name("dirichlettimedependent");
} // constructor


// ----------------------------------------------------------------------
// Destructor.
pylith::bc::DirichletTimeDependent::~DirichletTimeDependent(void)
{ // destructor
    deallocate();
} // destructor


// ----------------------------------------------------------------------
// Deallocate PETSc and local data structures.
void
pylith::bc::DirichletTimeDependent::deallocate(void)
{ // deallocate
    PYLITH_METHOD_BEGIN;

    DirichletNew::deallocate();

    PYLITH_METHOD_END;
} // deallocate


// ----------------------------------------------------------------------
// Use initial value term in time history expression.
void
pylith::bc::DirichletTimeDependent::useInitial(const bool value)
{ // useInitial
    PYLITH_JOURNAL_DEBUG("useInitial(value="<<value<<")");

    _useInitial = value;
} // useInitial


// ----------------------------------------------------------------------
// Get flag associated with using initial value term in time history expression.
bool
pylith::bc::DirichletTimeDependent::useInitial(void) const
{ // useInitial
    return _useInitial;
} // useInitial


// ----------------------------------------------------------------------
// Use rate value term in time history expression.
void
pylith::bc::DirichletTimeDependent::useRate(const bool value)
{ // useRate
    PYLITH_JOURNAL_DEBUG("useRate(value="<<value<<")");

    _useRate = value;
} // useRate


// ----------------------------------------------------------------------
// Get flag associated with using rate value term in time history expression.
bool
pylith::bc::DirichletTimeDependent::useRate(void) const
{ // useRate
    return _useRate;
} // useRate


// ----------------------------------------------------------------------
// Use time history term in time history expression.
void
pylith::bc::DirichletTimeDependent::useTimeHistory(const bool value)
{ // useTimeHistory
    PYLITH_JOURNAL_DEBUG("useTimeHistory(value="<<value<<")");

    _useTimeHistory = value;
} // useTimeHistory


// ----------------------------------------------------------------------
// Get flag associated with using time history term in time history expression.
bool
pylith::bc::DirichletTimeDependent::useTimeHistory(void) const
{ // useTimeHistory
    return _useTimeHistory;
} // useTimeHistory


// ----------------------------------------------------------------------
// Setup auxiliary subfields (discretization and query fns).
void
pylith::bc::DirichletTimeDependent::_auxFieldsSetup(void)
{ // _auxFieldsSetup
    PYLITH_JOURNAL_DEBUG("_auxFieldsSetup()");

    PYLITH_JOURNAL_ERROR(":TODO: @brad Implement _auxFieldsSetup().");
}     // _auxFieldsSetup


// ----------------------------------------------------------------------
// Set kernels for RHS residual G(t,s).
void
pylith::bc::DirichletTimeDependent::_setFEKernelsConstraint(const topology::Field& solution)
{ // _setFEKernelsConstraint
    PYLITH_JOURNAL_DEBUG("setFEKernelsConstraint(solution="<<solution.label()<<")");

    PYLITH_JOURNAL_ERROR(":TODO: @brad Implement _setFEKernelsConstraint().");
} // _setFEKernelsConstraint

// End of file
