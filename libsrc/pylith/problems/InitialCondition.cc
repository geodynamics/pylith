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
// Copyright (c) 2010-2015 University of California, Davis
//
// See COPYING for license information.
//
// ======================================================================
//

#include <portinfo>

#include "InitialCondition.hh" // implementation of class methods

#include "pylith/topology/Field.hh" // USES Field

#include "pylith/utils/error.hh" // USES PYLITH_CHECK_ERROR
#include "pylith/utils/journals.hh" // USES PYLITH_COMPONENT_

// ----------------------------------------------------------------------
// Constructor
pylith::problems::InitialCondition::InitialCondition(void)
{}


// ---------------------------------------------------------------------------------------------------------------------
// Destructor
pylith::problems::InitialCondition::~InitialCondition(void) {
    deallocate();
} // destructor


// ---------------------------------------------------------------------------------------------------------------------
// Deallocate PETSc and local data structures.
void
pylith::problems::InitialCondition::deallocate(void) {
    PYLITH_METHOD_BEGIN;

    PYLITH_METHOD_END;
} // deallocate


// ---------------------------------------------------------------------------------------------------------------------
// Set fields for initial condition.
void
pylith::problems::InitialCondition::setFields(const char* fields[],
                                              const int numFields) {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("setFields(fields="<<fields<<", numFields="<<numFields<<")");

    if (numFields > 0) {
        _fields.resize(numFields);
        for (int i = 0; i < numFields; ++i) {
            _fields[i] = fields[i];
        } // for
    } // if

    PYLITH_METHOD_END;
} // setFields


// ---------------------------------------------------------------------------------------------------------------------
// Verify configuration is acceptable.
void
pylith::problems::InitialCondition::verifyConfiguration(const pylith::topology::Field& solution) const {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("verifyConfiguration(solution="<<solution.label()<<") empty method");

    PYLITH_METHOD_END;
} // verityConfiguration


// End of file
