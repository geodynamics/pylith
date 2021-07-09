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
// Copyright (c) 2010-2021 University of California, Davis
//
// See LICENSE.md for license information.
//
// ----------------------------------------------------------------------
//

#include <portinfo>

#include "DerivedFactoryElasticity.hh" // implementation of object methods

#include "pylith/topology/Field.hh" // USES Field

#include "spatialdata/units/Nondimensional.hh" // USES Nondimensional

#include "pylith/utils/error.hh" // USES PYLITH_METHOD*
#include "pylith/utils/journals.hh" // USES PYLITH_JOURNAL*

#include <cassert>

// ---------------------------------------------------------------------------------------------------------------------
// Default constructor.
pylith::materials::DerivedFactoryElasticity::DerivedFactoryElasticity(void) {
    GenericComponent::setName("derivedfactoryelasticity");
} // constructor


// ---------------------------------------------------------------------------------------------------------------------
// Destructor.
pylith::materials::DerivedFactoryElasticity::~DerivedFactoryElasticity(void) {}


// ---------------------------------------------------------------------------------------------------------------------
// Add Cauchy stress subfield to derived field.
void
pylith::materials::DerivedFactoryElasticity::addCauchyStress(void) {
    PYLITH_METHOD_BEGIN;
    PYLITH_JOURNAL_DEBUG("addCauchyStress(void)");

    const char* fieldName = "cauchy_stress";
    const char* componentNames[6] = { "stress_xx", "stress_yy", "stress_zz", "stress_xy", "stress_yz", "stress_xz" };
    const int stressSize = (3 == _spaceDim) ? 6 : (2 == _spaceDim) ? 4 : 1;
    const PylithReal pressureScale = _normalizer->getPressureScale();

    pylith::topology::Field::Description description;
    description.label = fieldName;
    description.alias = fieldName;
    description.vectorFieldType = (3 == _spaceDim) ? pylith::topology::Field::TENSOR : pylith::topology::Field::OTHER;
    description.numComponents = stressSize;
    description.componentNames.resize(stressSize);
    for (int i = 0; i < stressSize; ++i) {
        description.componentNames[i] = componentNames[i];
    } // for
    description.scale = pressureScale;
    description.validator = NULL;

    _field->subfieldAdd(description, getSubfieldDiscretization(fieldName));

    PYLITH_METHOD_END;
} // addCauchyStress


// ---------------------------------------------------------------------------------------------------------------------
// Add Cauchy (infinitesimal) strain subfield to derived fields.
void
pylith::materials::DerivedFactoryElasticity::addCauchyStrain(void) {
    PYLITH_METHOD_BEGIN;
    PYLITH_JOURNAL_DEBUG("addCauchyStrain(void)");

    const char* fieldName = "cauchy_strain";
    const char* componentNames[6] = { "strain_xx", "strain_yy", "strain_zz", "strain_xy", "strain_yz", "strain_xz" };
    const int strainSize = (3 == _spaceDim) ? 6 : (2 == _spaceDim) ? 4 : 1;

    pylith::topology::Field::Description description;
    description.label = fieldName;
    description.alias = fieldName;
    description.vectorFieldType = (3 == _spaceDim) ? pylith::topology::Field::TENSOR : pylith::topology::Field::OTHER;
    description.numComponents = strainSize;
    description.componentNames.resize(strainSize);
    for (int i = 0; i < strainSize; ++i) {
        description.componentNames[i] = componentNames[i];
    } // for
    description.scale = 1.0;
    description.validator = NULL;

    _field->subfieldAdd(description, getSubfieldDiscretization(fieldName));

    PYLITH_METHOD_END;
} // addCauchyStrain


// ---------------------------------------------------------------------------------------------------------------------
// Add subfields using discretizations provided.
void
pylith::materials::DerivedFactoryElasticity::addSubfields(void) {
    PYLITH_METHOD_BEGIN;
    PYLITH_JOURNAL_DEBUG("addSubfields(void)");

    if (_subfieldDiscretizations.find("cauchy_stress") != _subfieldDiscretizations.end()) {
        addCauchyStress();
    } // if
    if (_subfieldDiscretizations.find("cauchy_strain") != _subfieldDiscretizations.end()) {
        addCauchyStrain();
    } // if

    PYLITH_METHOD_END;

}


// End of file
