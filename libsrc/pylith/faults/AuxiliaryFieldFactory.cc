// =================================================================================================
// This code is part of PyLith, developed through the Computational Infrastructure
// for Geodynamics (https://github.com/geodynamics/pylith).
//
// Copyright (c) 2010-2023, University of California, Davis and the PyLith Development Team.
// All rights reserved.
//
// See https://mit-license.org/ and LICENSE.md and for license information.
// =================================================================================================

#include <portinfo>

#include "pylith/faults/AuxiliaryFieldFactory.hh" // implementation of object methods

#include "pylith/topology/Field.hh" // USES Field
#include "pylith/topology/FieldQuery.hh" // HOLDSA FieldQuery

#include "spatialdata/units/Nondimensional.hh" // USES Nondimensional

#include "pylith/utils/error.hh" // USES PYLITH_METHOD*
#include "pylith/utils/journals.hh" // USES PYLITH_JOURNAL*

#include <cmath> // USES pow()
#include <cassert>

// ------------------------------------------------------------------------------------------------
// Default constructor.
pylith::faults::AuxiliaryFieldFactory::AuxiliaryFieldFactory(void) {
    GenericComponent::setName("auxiliaryfactorykinematic");
} // constructor


// ------------------------------------------------------------------------------------------------
// Destructor.
pylith::faults::AuxiliaryFieldFactory::~AuxiliaryFieldFactory(void) {}


// ------------------------------------------------------------------------------------------------
// Add fault slip subfield to auxiliary fields.
void
pylith::faults::AuxiliaryFieldFactory::addSlip(void) {
    PYLITH_METHOD_BEGIN;
    PYLITH_JOURNAL_DEBUG("addSlip(void)");

    const char* fieldName = "slip";
    const char* componentNames[3] = { "slip_opening", "slip_left_lateral", "slip_reverse" };

    const PylithReal lengthScale = _normalizer->getLengthScale();

    pylith::topology::Field::Description description;
    description.label = fieldName;
    description.alias = fieldName;
    description.vectorFieldType = pylith::topology::Field::VECTOR;
    description.numComponents = _spaceDim;
    description.componentNames.resize(_spaceDim);
    for (int i = 0; i < _spaceDim; ++i) {
        description.componentNames[i] = componentNames[i];
    } // for
    description.scale = lengthScale;
    description.validator = NULL;

    _field->subfieldAdd(description, getSubfieldDiscretization(fieldName));
    // No query; populated by kinematic source at beginning of time step.

    PYLITH_METHOD_END;
} // addSlip


// ------------------------------------------------------------------------------------------------
// Add fault slip rate subfield to auxiliary fields.
void
pylith::faults::AuxiliaryFieldFactory::addSlipRate(void) {
    PYLITH_METHOD_BEGIN;
    PYLITH_JOURNAL_DEBUG("addSlipRate(void)");

    const char* fieldName = "slip_rate";
    const char* componentNames[3] = { "slip_rate_opening", "slip_rate_left_lateral", "slip_rate_reverse" };

    const PylithReal velocityScale = _normalizer->getLengthScale() / _normalizer->getTimeScale();

    pylith::topology::Field::Description description;
    description.label = fieldName;
    description.alias = fieldName;
    description.vectorFieldType = pylith::topology::Field::VECTOR;
    description.numComponents = _spaceDim;
    description.componentNames.resize(_spaceDim);
    for (int i = 0; i < _spaceDim; ++i) {
        description.componentNames[i] = componentNames[i];
    } // for
    description.scale = velocityScale;
    description.validator = NULL;

    _field->subfieldAdd(description, getSubfieldDiscretization(fieldName));
    // No query; populated by kinematic source at beginning of time step.

    PYLITH_METHOD_END;
} // addSlipRate


// ------------------------------------------------------------------------------------------------
// Add fault slip acceleration subfield to auxiliary fields.
void
pylith::faults::AuxiliaryFieldFactory::addSlipAcceleration(void) {
    PYLITH_METHOD_BEGIN;
    PYLITH_JOURNAL_DEBUG("addSlipAcceleration(void)");

    const char* fieldName = "slip_acceleration";
    const char* componentNames[3] = {
        "slip_acceleration_opening",
        "slip_acceleration_left_lateral",
        "slip_acceleration_reverse",
    };

    const PylithReal accelerationScale = _normalizer->getLengthScale() / pow(_normalizer->getTimeScale(), 2);

    pylith::topology::Field::Description description;
    description.label = fieldName;
    description.alias = fieldName;
    description.vectorFieldType = pylith::topology::Field::VECTOR;
    description.numComponents = _spaceDim;
    description.componentNames.resize(_spaceDim);
    for (int i = 0; i < _spaceDim; ++i) {
        description.componentNames[i] = componentNames[i];
    } // for
    description.scale = accelerationScale;
    description.validator = NULL;

    _field->subfieldAdd(description, getSubfieldDiscretization(fieldName));
    // No query; populated by kinematic source at beginning of time step.

    PYLITH_METHOD_END;
} // addSlipAcc


// End of file
