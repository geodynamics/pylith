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
pylith::problems::InitialCondition::setSubfields(const char* subfields[],
                                                 const int numSubfields) {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("setFields(subfields="<<subfields<<", numSubfields="<<numSubfields<<")");

    if (numSubfields > 0) {
        _subfields.resize(numSubfields);
        for (int i = 0; i < numSubfields; ++i) {
            _subfields[i] = subfields[i];
        } // for
    } // if

    PYLITH_METHOD_END;
} // setFields


// ---------------------------------------------------------------------------------------------------------------------
// Verify configuration is acceptable.
void
pylith::problems::InitialCondition::verifyConfiguration(const pylith::topology::Field& solution) const {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("verifyConfiguration(solution="<<solution.getLabel()<<")");

    const size_t numSubfields = _subfields.size();
    for (size_t i = 0; i < numSubfields; ++i) {
        if (!solution.hasSubfield(_subfields[i].c_str())) {
            std::ostringstream msg;
            msg << "Cannot specify initial conditions for solution subfield '"<< _subfields[i]
                << "' in component '" << PyreComponent::getIdentifier() << "'"
                << "; field is not in solution.";
            throw std::runtime_error(msg.str());
        } // if
    } // for

    PYLITH_METHOD_END;
} // verityConfiguration


// End of file
