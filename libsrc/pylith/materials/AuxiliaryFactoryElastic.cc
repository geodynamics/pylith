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
// Copyright (c) 2010-2015 University of California, Davis
//
// See COPYING for license information.
//
// ----------------------------------------------------------------------
//

#include <portinfo>

#include "AuxiliaryFactoryElastic.hh" // implementation of object methods

#include "Material.hh" // USES Material
#include "Query.hh" // USES Query

#include "pylith/topology/Field.hh" // USES Field
#include "pylith/topology/FieldQuery.hh" // HOLDSA FieldQuery

#include "spatialdata/units/Nondimensional.hh" // USES Nondimensional
#include "spatialdata/spatialdb/GravityField.hh" // USES GravityField

#include "pylith/utils/journals.hh" // USES PYLITH_JOURNAL*

#include <cassert>

// ---------------------------------------------------------------------------------------------------------------------
// Default constructor.
pylith::materials::AuxiliaryFactoryElastic::AuxiliaryFactoryElastic(void) {
    GenericComponent::name("auxiliaryfactoryelastic");
} // constructor


// ---------------------------------------------------------------------------------------------------------------------
// Destructor.
pylith::materials::AuxiliaryFactoryElastic::~AuxiliaryFactoryElastic(void) {}


// ---------------------------------------------------------------------------------------------------------------------
// Add density subfield to auxiliary fields.
void
pylith::materials::AuxiliaryFactoryElastic::addDensity(void) {
    PYLITH_METHOD_BEGIN;
    PYLITH_JOURNAL_DEBUG("addDensity(void)");

    const char* fieldName = "density";
    const PylithReal densityScale = _normalizer->densityScale();

    pylith::topology::Field::Description description;
    description.label = fieldName;
    description.alias = fieldName;
    description.vectorFieldType = pylith::topology::Field::SCALAR;
    description.numComponents = 1;
    description.componentNames.resize(1);
    description.componentNames[0] = fieldName;
    description.scale = densityScale;
    description.validator = pylith::topology::FieldQuery::validatorPositive;

    _field->subfieldAdd(description, getSubfieldDiscretization(fieldName));
    _setSubfieldQueryFn(fieldName, pylith::topology::FieldQuery::dbQueryGeneric);

    PYLITH_METHOD_END;
} // addDensity


// ---------------------------------------------------------------------------------------------------------------------
// Add shear modulus subfield to auxiliary fields.
void
pylith::materials::AuxiliaryFactoryElastic::addShearModulus(void) {
    PYLITH_METHOD_BEGIN;
    PYLITH_JOURNAL_DEBUG("addShearModulus(void)");

    const char* fieldName = "shear_modulus";
    const PylithReal pressureScale = _normalizer->pressureScale();

    pylith::topology::Field::Description description;
    description.label = fieldName;
    description.alias = fieldName;
    description.vectorFieldType = pylith::topology::Field::SCALAR;
    description.numComponents = 1;
    description.componentNames.resize(1);
    description.componentNames[0] = fieldName;
    description.scale = pressureScale;
    description.validator = NULL;

    _field->subfieldAdd(description, getSubfieldDiscretization(fieldName));
    _setSubfieldQueryFn(fieldName, pylith::materials::Query::dbQueryShearModulus);

    PYLITH_METHOD_END;
} // addShearModulus


// ---------------------------------------------------------------------------------------------------------------------
// Add bulk modulus subfield to auxiliary fields.
void
pylith::materials::AuxiliaryFactoryElastic::addBulkModulus(void) {
    PYLITH_METHOD_BEGIN;
    PYLITH_JOURNAL_DEBUG("addBulkModulus(void)");

    const char* fieldName = "bulk_modulus";
    const PylithReal pressureScale = _normalizer->pressureScale();

    pylith::topology::Field::Description description;
    description.label = fieldName;
    description.alias = fieldName;
    description.vectorFieldType = pylith::topology::Field::SCALAR;
    description.numComponents = 1;
    description.componentNames.resize(1);
    description.componentNames[0] = fieldName;
    description.scale = pressureScale;
    description.validator = NULL;

    _field->subfieldAdd(description, getSubfieldDiscretization(fieldName));
    _setSubfieldQueryFn(fieldName, pylith::materials::Query::dbQueryBulkModulus);

    PYLITH_METHOD_END;
} // addBulkModulus


// ---------------------------------------------------------------------------------------------------------------------
// Add gravity subfield to auxiliary fields.
void
pylith::materials::AuxiliaryFactoryElastic::addGravityField(spatialdata::spatialdb::GravityField* gf) {
    PYLITH_METHOD_BEGIN;
    PYLITH_JOURNAL_DEBUG("addGravityField(void)");

    const char* fieldName = "gravitational_acceleration";
    const char* componentNames[3] = { "gravitational_acceleration_x", "gravitational_acceleration_y", "gravitational_acceleration_z" };

    const PylithReal lengthScale = _normalizer->lengthScale();
    const PylithReal timeScale = _normalizer->timeScale();
    const PylithReal accelerationScale = lengthScale / (timeScale * timeScale);

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
    _setSubfieldQueryFn(fieldName, pylith::materials::Query::dbQueryGravityField, gf);

    PYLITH_METHOD_END;
} // addGravityField


// ---------------------------------------------------------------------------------------------------------------------
// Add body force subfield to auxiliary fields.
void
pylith::materials::AuxiliaryFactoryElastic::addBodyForce(void) {
    PYLITH_METHOD_BEGIN;
    PYLITH_JOURNAL_DEBUG("addBodyForce(void)");

    const char* fieldName = "body_force";
    const char* componentNames[3] = { "body_force_x", "body_force_y", "body_force_z" };

    const PylithReal forceScale = _normalizer->pressureScale() / _normalizer->lengthScale();

    pylith::topology::Field::Description description;
    description.label = fieldName;
    description.alias = fieldName;
    description.vectorFieldType = pylith::topology::Field::VECTOR;
    description.numComponents = _spaceDim;
    description.componentNames.resize(_spaceDim);
    for (int i = 0; i < _spaceDim; ++i) {
        description.componentNames[i] = componentNames[i];
    } // for
    description.scale = forceScale;
    description.validator = NULL;

    _field->subfieldAdd(description, getSubfieldDiscretization(fieldName));
    _setSubfieldQueryFn(fieldName, pylith::topology::FieldQuery::dbQueryGeneric);

    PYLITH_METHOD_END;
} // addBodyForce


// ---------------------------------------------------------------------------------------------------------------------
// Add reference stress subfield to auxiliary fields.
void
pylith::materials::AuxiliaryFactoryElastic::addReferenceStress(void) {
    PYLITH_METHOD_BEGIN;
    PYLITH_JOURNAL_DEBUG("addReferenceStress(void)");

    const char* fieldName = "reference_stress";
    const char* componentNames[6] = { "reference_stress_xx", "reference_stress_yy", "reference_stress_zz", "reference_stress_xy", "reference_stress_yz", "reference_stress_xz" };
    const int stressSize = (3 == _spaceDim) ? 6 : (2 == _spaceDim) ? 4 : 1;
    const PylithReal pressureScale = _normalizer->pressureScale();

    pylith::topology::Field::Description description;
    description.label = fieldName;
    description.alias = fieldName;
    description.vectorFieldType = pylith::topology::Field::OTHER;
    description.numComponents = stressSize;
    description.componentNames.resize(stressSize);
    for (int i = 0; i < stressSize; ++i) {
        description.componentNames[i] = componentNames[i];
    } // for
    description.scale = pressureScale;
    description.validator = NULL;

    _field->subfieldAdd(description, getSubfieldDiscretization(fieldName));
    _setSubfieldQueryFn(fieldName, pylith::topology::FieldQuery::dbQueryGeneric);

    PYLITH_METHOD_END;
} // addReferenceStress


// ---------------------------------------------------------------------------------------------------------------------
// Add reference strain subfield to auxiliary fields.
void
pylith::materials::AuxiliaryFactoryElastic::addReferenceStrain(void) {
    PYLITH_METHOD_BEGIN;
    PYLITH_JOURNAL_DEBUG("addRefrenceStrain(void)");

    const char* fieldName = "reference_strain";
    const char* componentNames[6] = { "reference_strain_xx", "reference_strain_yy", "reference_strain_zz", "reference_strain_xy", "reference_strain_yz", "reference_strain_xz" };
    const int strainSize = (3 == _spaceDim) ? 6 : (2 == _spaceDim) ? 4 : 1;

    pylith::topology::Field::Description description;
    description.label = fieldName;
    description.alias = fieldName;
    description.vectorFieldType = pylith::topology::Field::OTHER;
    description.numComponents = strainSize;
    description.componentNames.resize(strainSize);
    for (int i = 0; i < strainSize; ++i) {
        description.componentNames[i] = componentNames[i];
    } // for
    description.scale = 1.0;
    description.validator = NULL;

    _field->subfieldAdd(description, getSubfieldDiscretization(fieldName));
    _setSubfieldQueryFn(fieldName, pylith::topology::FieldQuery::dbQueryGeneric);

    PYLITH_METHOD_END;
} // addReferenceStrain


// End of file
