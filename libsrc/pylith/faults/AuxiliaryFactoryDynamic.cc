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
// Copyright (c) 2010-2022 University of California, Davis
//
// See LICENSE.md.md for license information.
//
// ----------------------------------------------------------------------
//

#include <portinfo>

#include "AuxiliaryFactoryDynamic.hh" // implementation of object methods

#include "pylith/topology/Field.hh" // USES Field
#include "pylith/topology/FieldQuery.hh" // HOLDSA FieldQuery

#include "spatialdata/units/Nondimensional.hh" // USES Nondimensional

#include "pylith/utils/error.hh" // USES PYLITH_METHOD*
#include "pylith/utils/journals.hh" // USES PYLITH_JOURNAL*

#include <cmath> // USES pow()
#include <cassert>

// ------------------------------------------------------------------------------------------------
// Default constructor.
pylith::faults::AuxiliaryFactoryDynamic::AuxiliaryFactoryDynamic(void) {
    GenericComponent::setName("auxiliaryfactorydynamic");
} // constructor


// ------------------------------------------------------------------------------------------------
// Destructor.
pylith::faults::AuxiliaryFactoryDynamic::~AuxiliaryFactoryDynamic(void) {}


// ------------------------------------------------------------------------------------------------
// Add fault traction perturbation to auxiliary fields.
void
pylith::faults::AuxiliaryFactoryDynamic::addTractionPerturbation(void) {
    PYLITH_METHOD_BEGIN;
    PYLITH_JOURNAL_DEBUG("addTractionPerturbation(void)");

    const char* fieldName = "perturbation_traction";
    const char* componentNames[3] = { "traction_opening", "traction_left_lateral", "traction_reverse" };

    const PylithReal pressureScale = _normalizer->getPressureScale();

    pylith::topology::Field::Description description;
    description.label = fieldName;
    description.alias = fieldName;
    description.vectorFieldType = pylith::topology::Field::VECTOR;
    description.numComponents = _spaceDim;
    description.componentNames.resize(_spaceDim);
    for (int i = 0; i < _spaceDim; ++i) {
        description.componentNames[i] = componentNames[i];
    } // for
    description.scale = pressureScale;
    description.validator = NULL;

    _field->subfieldAdd(description, getSubfieldDiscretization(fieldName));
    // No query; populated by traction perturbation at beginning of time step.

    PYLITH_METHOD_END;
} // addTractionPerturbation


// ------------------------------------------------------------------------------------------------
// Add fault rheology traction subfield to auxiliary fields.
void
pylith::faults::AuxiliaryFactoryDynamic::addTractionRheology(void) {
    PYLITH_METHOD_BEGIN;
    PYLITH_JOURNAL_DEBUG("addTractionRheology(void)");

    const char* fieldName = "rheology_traction";
    const char* componentNames[3] = { "traction_opening", "traction_left_lateral", "traction_reverse" };

    const PylithReal pressureScale = _normalizer->getPressureScale();

    pylith::topology::Field::Description description;
    description.label = fieldName;
    description.alias = fieldName;
    description.vectorFieldType = pylith::topology::Field::VECTOR;
    description.numComponents = _spaceDim;
    description.componentNames.resize(_spaceDim);
    for (int i = 0; i < _spaceDim; ++i) {
        description.componentNames[i] = componentNames[i];
    } // for
    description.scale = pressureScale;
    description.validator = NULL;

    _field->subfieldAdd(description, getSubfieldDiscretization(fieldName));
    // No query; populated by dynamic source at beginning of time step.

    PYLITH_METHOD_END;
} // addTractionRheology


// End of file
