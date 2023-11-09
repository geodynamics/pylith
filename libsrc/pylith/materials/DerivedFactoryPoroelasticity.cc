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

#include "DerivedFactoryPoroelasticity.hh" // implementation of object methods

#include "pylith/topology/Field.hh" // USES Field

#include "pylith/topology/FieldQuery.hh" // HOLDSA FieldQuery

#include "spatialdata/units/Nondimensional.hh" // USES Nondimensional

#include "pylith/utils/error.hh" // USES PYLITH_METHOD*
#include "pylith/utils/journals.hh" // USES PYLITH_JOURNAL*

#include <cassert>

// ------------------------------------------------------------------------------------------------
// Default constructor.
pylith::materials::DerivedFactoryPoroelasticity::DerivedFactoryPoroelasticity(void) {
    GenericComponent::setName("derivedfactoryporoelasticity");
} // constructor


// ------------------------------------------------------------------------------------------------
// Destructor.
pylith::materials::DerivedFactoryPoroelasticity::~DerivedFactoryPoroelasticity(void) {}


// ------------------------------------------------------------------------------------------------
// Add bulk density subfield to derived fields
void
pylith::materials::DerivedFactoryPoroelasticity::addBulkDensity(void) {
    PYLITH_METHOD_BEGIN;
    PYLITH_JOURNAL_DEBUG("addBulkDensity(void)");

    const char* fieldName = "bulk_density";
    const PylithReal densityScale = _normalizer->getDensityScale();

    
    pylith::topology::Field::Description description;
    description.label = fieldName;
    description.alias = fieldName;
    description.vectorFieldType = pylith::topology::Field::SCALAR;
    description.numComponents = 1;
    description.componentNames.resize(1);
    description.componentNames[0] = fieldName;
    description.scale = densityScale;
    description.validator = pylith::topology::FieldQuery::validatorNonnegative;;

    _field->subfieldAdd(description, getSubfieldDiscretization(fieldName));

    PYLITH_METHOD_END;
} //addBulkDensity

// ------------------------------------------------------------------------------------------------
// Add water content subfield to derived fields
void
pylith::materials::DerivedFactoryPoroelasticity::addWaterContent(void) {
    PYLITH_METHOD_BEGIN;
    PYLITH_JOURNAL_DEBUG("addWaterContent(void)");

    const char* fieldName = "water_content";
    
    pylith::topology::Field::Description description;
    description.label = fieldName;
    description.alias = fieldName;
    description.vectorFieldType = pylith::topology::Field::SCALAR;
    description.numComponents = 1;
    description.componentNames.resize(1);
    description.componentNames[0] = fieldName;
    description.scale = 1.0;

    _field->subfieldAdd(description, getSubfieldDiscretization(fieldName));

    PYLITH_METHOD_END;
} //addWaterContent

// ------------------------------------------------------------------------------------------------
// Add subfields using discretizations provided.
void
pylith::materials::DerivedFactoryPoroelasticity::addSubfields(void) {
    PYLITH_METHOD_BEGIN;
    PYLITH_JOURNAL_DEBUG("addSubfields(void)");

    if (_subfieldDiscretizations.find("cauchy_stress") != _subfieldDiscretizations.end()) {
        addCauchyStress();
    } // if
    if (_subfieldDiscretizations.find("cauchy_strain") != _subfieldDiscretizations.end()) {
        addCauchyStrain();
    } // if
    if (_subfieldDiscretizations.find("bulk_density") != _subfieldDiscretizations.end()) {
        addBulkDensity();
    } // if
    if (_subfieldDiscretizations.find("water_content") != _subfieldDiscretizations.end()) {
        addWaterContent();
    } // if

    PYLITH_METHOD_END;
} //addSubfields

// End of file
