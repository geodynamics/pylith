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
// See LICENSE.md for license information.
//
// ----------------------------------------------------------------------
//

#include <portinfo>

#include "DerivedFactoryKinematic.hh" // implementation of object methods

#include "pylith/topology/Field.hh" // USES Field

#include "spatialdata/units/Nondimensional.hh" // USES Nondimensional

#include "pylith/utils/error.hh" // USES PYLITH_METHOD*
#include "pylith/utils/journals.hh" // USES PYLITH_JOURNAL*

#include <cassert>

// ------------------------------------------------------------------------------------------------
// Default constructor.
pylith::faults::DerivedFactoryKinematic::DerivedFactoryKinematic(void) {
    GenericComponent::setName("derivedfactorykinematic");
} // constructor


// ------------------------------------------------------------------------------------------------
// Destructor.
pylith::faults::DerivedFactoryKinematic::~DerivedFactoryKinematic(void) {}


// ------------------------------------------------------------------------------------------------
// Add Cauchy stress subfield to derived field.
void
pylith::faults::DerivedFactoryKinematic::addTractionChange(void) {
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
pylith::faults::DerivedFactoryKinematic::addSubfields(void) {
    PYLITH_METHOD_BEGIN;
    PYLITH_JOURNAL_DEBUG("addSubfields(void)");

    if (_subfieldDiscretizations.find("traction_change") != _subfieldDiscretizations.end()) {
        addTractionChange();
    } // if

    PYLITH_METHOD_END;

}


// End of file
