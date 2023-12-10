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

#include "pylith/faults/DerivedFieldFactory.hh" // implementation of object methods

#include "pylith/topology/Field.hh" // USES Field

#include "spatialdata/units/Nondimensional.hh" // USES Nondimensional

#include "pylith/utils/error.hh" // USES PYLITH_METHOD*
#include "pylith/utils/journals.hh" // USES PYLITH_JOURNAL*

#include <cassert>

// ------------------------------------------------------------------------------------------------
// Default constructor.
pylith::faults::DerivedFieldFactory::DerivedFieldFactory(void) {
    GenericComponent::setName("derivedfactorykinematic");
} // constructor


// ------------------------------------------------------------------------------------------------
// Destructor.
pylith::faults::DerivedFieldFactory::~DerivedFieldFactory(void) {}


// ------------------------------------------------------------------------------------------------
// Add Cauchy stress subfield to derived field.
void
pylith::faults::DerivedFieldFactory::addTractionChange(void) {
    PYLITH_METHOD_BEGIN;
    PYLITH_JOURNAL_DEBUG("addTractionChange(void)");

    const char* fieldName = "traction_change";
    const char* componentNames[3] = { "traction_change_opening", "traction_change_left_lateral", "traction_change_reverse" };
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

    PYLITH_METHOD_END;
} // addCauchyStress


// ------------------------------------------------------------------------------------------------
// Add subfields using discretizations provided.
void
pylith::faults::DerivedFieldFactory::addSubfields(void) {
    PYLITH_METHOD_BEGIN;
    PYLITH_JOURNAL_DEBUG("addSubfields(void)");

    if (_subfieldDiscretizations.find("traction_change") != _subfieldDiscretizations.end()) {
        addTractionChange();
    } // if

    PYLITH_METHOD_END;

}


// End of file
