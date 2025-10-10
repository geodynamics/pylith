// =================================================================================================
// This code is part of PyLith, developed through the Computational Infrastructure
// for Geodynamics (https://github.com/geodynamics/pylith).
//
// Copyright (c) 2010-2025, University of California, Davis and the PyLith Development Team.
// All rights reserved.
//
// See https://mit-license.org/ and LICENSE.md and for license information.
// =================================================================================================

#include <portinfo>

#include "pylith/materials/AuxiliaryFactoryElastic.hh" // implementation of object methods

#include "pylith/materials/Query.hh" // USES Query

#include "pylith/topology/Field.hh" // USES Field
#include "pylith/topology/FieldQuery.hh" // HOLDSA FieldQuery

#include "pylith/scales/Scales.hh" // USES Scales
#include "pylith/scales/ElasticityScales.hh" // USES ElasticityScales

#include "pylith/utils/error.hh" // USES PYLITH_METHOD*
#include "pylith/utils/journals.hh" // USES PYLITH_JOURNAL*

#include <cassert>

// ------------------------------------------------------------------------------------------------
// Default constructor.
pylith::materials::AuxiliaryFactoryElastic::AuxiliaryFactoryElastic(void) {
    GenericComponent::setName("auxiliaryfactoryelastic");
} // constructor


// ------------------------------------------------------------------------------------------------
// Destructor.
pylith::materials::AuxiliaryFactoryElastic::~AuxiliaryFactoryElastic(void) {}


// ------------------------------------------------------------------------------------------------
// Add shear modulus subfield to auxiliary fields.
void
pylith::materials::AuxiliaryFactoryElastic::addShearModulus(void) {
    PYLITH_METHOD_BEGIN;
    PYLITH_JOURNAL_DEBUG("addShearModulus(void)");

    const char* subfieldName = "shear_modulus";
    const PylithReal rigidityScale = _scales->getRigidityScale();

    pylith::topology::Field::Description description;
    description.label = subfieldName;
    description.alias = subfieldName;
    description.vectorFieldType = pylith::topology::Field::SCALAR;
    description.numComponents = 1;
    description.componentNames.resize(1);
    description.componentNames[0] = subfieldName;
    description.scale = rigidityScale;
    description.validator = pylith::topology::FieldQuery::validatorNonnegative;

    _field->subfieldAdd(description, getSubfieldDiscretization(subfieldName));
    pylith::materials::Query::shearModulusFromVM(subfieldName, this);

    PYLITH_METHOD_END;
} // addShearModulus


// ------------------------------------------------------------------------------------------------
// Add bulk modulus subfield to auxiliary fields.
void
pylith::materials::AuxiliaryFactoryElastic::addBulkModulus(void) {
    PYLITH_METHOD_BEGIN;
    PYLITH_JOURNAL_DEBUG("addBulkModulus(void)");

    const char* subfieldName = "bulk_modulus";
    const PylithReal rigidityScale = _scales->getRigidityScale();

    pylith::topology::Field::Description description;
    description.label = subfieldName;
    description.alias = subfieldName;
    description.vectorFieldType = pylith::topology::Field::SCALAR;
    description.numComponents = 1;
    description.componentNames.resize(1);
    description.componentNames[0] = subfieldName;
    description.scale = rigidityScale;
    description.validator = pylith::topology::FieldQuery::validatorPositive;

    _field->subfieldAdd(description, getSubfieldDiscretization(subfieldName));
    pylith::materials::Query::bulkModulusFromVM(subfieldName, this);

    PYLITH_METHOD_END;
} // addBulkModulus


// ------------------------------------------------------------------------------------------------
// Add reference stress subfield to auxiliary fields.
void
pylith::materials::AuxiliaryFactoryElastic::addReferenceStress(void) {
    PYLITH_METHOD_BEGIN;
    PYLITH_JOURNAL_DEBUG("addReferenceStress(void)");

    const char* subfieldName = "reference_stress";
    const char* componentNames[6] = {
        "reference_stress_xx",
        "reference_stress_yy",
        "reference_stress_zz",
        "reference_stress_xy",
        "reference_stress_yz",
        "reference_stress_xz" };
    const int stressSize = (3 == _spaceDim) ? 6 : (2 == _spaceDim) ? 4 : 1;
    const PylithReal stressScale = pylith::scales::ElasticityScales::getStressScale(*_scales);

    pylith::topology::Field::Description description;
    description.label = subfieldName;
    description.alias = subfieldName;
    description.vectorFieldType = (3 == _spaceDim) ? pylith::topology::Field::TENSOR : pylith::topology::Field::OTHER;
    description.numComponents = stressSize;
    description.componentNames.resize(stressSize);
    for (int i = 0; i < stressSize; ++i) {
        description.componentNames[i] = componentNames[i];
    } // for
    description.scale = stressScale;
    description.validator = NULL;

    _field->subfieldAdd(description, getSubfieldDiscretization(subfieldName));
    this->setSubfieldQuery(subfieldName);

    PYLITH_METHOD_END;
} // addReferenceStress


// ------------------------------------------------------------------------------------------------
// Add reference strain subfield to auxiliary fields.
void
pylith::materials::AuxiliaryFactoryElastic::addReferenceStrain(void) {
    PYLITH_METHOD_BEGIN;
    PYLITH_JOURNAL_DEBUG("addRefrenceStrain(void)");

    const char* subfieldName = "reference_strain";
    const char* componentNames[6] = {
        "reference_strain_xx",
        "reference_strain_yy",
        "reference_strain_zz",
        "reference_strain_xy",
        "reference_strain_yz",
        "reference_strain_xz"
    };
    const int strainSize = (3 == _spaceDim) ? 6 : (2 == _spaceDim) ? 4 : 1;
    const PylithReal strainScale = pylith::scales::ElasticityScales::getStrainScale(*_scales);

    pylith::topology::Field::Description description;
    description.label = subfieldName;
    description.alias = subfieldName;
    description.vectorFieldType = (3 == _spaceDim) ? pylith::topology::Field::TENSOR : pylith::topology::Field::OTHER;
    description.numComponents = strainSize;
    description.componentNames.resize(strainSize);
    for (int i = 0; i < strainSize; ++i) {
        description.componentNames[i] = componentNames[i];
    } // for
    description.scale = strainScale;
    description.validator = NULL;

    _field->subfieldAdd(description, getSubfieldDiscretization(subfieldName));
    this->setSubfieldQuery(subfieldName);

    PYLITH_METHOD_END;
} // addReferenceStrain


// End of file
