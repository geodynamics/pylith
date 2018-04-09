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

#include "AuxiliaryFactory.hh" // implementation of object methods

#include "Material.hh" // USES Material
#include "Query.hh" // USES Query

#include "pylith/topology/Field.hh" // USES Field
#include "pylith/topology/FieldQuery.hh" // HOLDSA FieldQuery

#include "spatialdata/units/Nondimensional.hh" // USES Nondimensional
#include "spatialdata/spatialdb/GravityField.hh" // USES GravityField

#include "pylith/utils/journals.hh" // USES PYLITH_JOURNAL*

#include <cassert>

// ----------------------------------------------------------------------
const char* pylith::materials::AuxiliaryFactory::_genericComponent = "materialauxiliaryfactory";

// ----------------------------------------------------------------------
// Default constructor.
pylith::materials::AuxiliaryFactory::AuxiliaryFactory(void)
{ // constructor
    GenericComponent::name(_genericComponent);
} // constructor

// ----------------------------------------------------------------------
// Destructor.
pylith::materials::AuxiliaryFactory::~AuxiliaryFactory(void)
{ // destructor
} // destructor

// ----------------------------------------------------------------------
// Add density subfield to auxiliary fields.
void
pylith::materials::AuxiliaryFactory::density(void)
{ // density
    PYLITH_METHOD_BEGIN;
    PYLITH_JOURNAL_DEBUG("density(void)");

    const char* fieldName = "density";
    const PylithReal densityScale = _normalizer->densityScale();

    pylith::topology::Field::Description description;
    description.label = fieldName;
    description.vectorFieldType = pylith::topology::Field::SCALAR;
    description.numComponents = 1;
    description.componentNames.resize(1);
    description.componentNames[0] = fieldName;
    description.scale = densityScale;
    description.validator = pylith::topology::FieldQuery::validatorPositive;

    _field->subfieldAdd(description, subfieldDiscretization(fieldName));
    _subfieldQueryFn(fieldName, pylith::topology::FieldQuery::dbQueryGeneric);

    PYLITH_METHOD_END;
} // density


// ----------------------------------------------------------------------
// Add shear modulus subfield to auxiliary fields.
void
pylith::materials::AuxiliaryFactory::shearModulus(void)
{ // shearModulus
    PYLITH_METHOD_BEGIN;
    PYLITH_JOURNAL_DEBUG("shearModulus(void)");

    const char* fieldName = "shear_modulus";
    const PylithReal pressureScale = _normalizer->pressureScale();

    pylith::topology::Field::Description description;
    description.label = fieldName;
    description.vectorFieldType = pylith::topology::Field::SCALAR;
    description.numComponents = 1;
    description.componentNames.resize(1);
    description.componentNames[0] = fieldName;
    description.scale = pressureScale;
    description.validator = NULL;

    _field->subfieldAdd(description, subfieldDiscretization(fieldName));
    _subfieldQueryFn(fieldName, pylith::materials::Query::dbQueryShearModulus);

    PYLITH_METHOD_END;
} // shearModulus


// ----------------------------------------------------------------------
// Add bulk modulus subfield to auxiliary fields.
void
pylith::materials::AuxiliaryFactory::bulkModulus(void)
{ // bulkModulus
    PYLITH_METHOD_BEGIN;
    PYLITH_JOURNAL_DEBUG("bulkModulus(void)");

    const char* fieldName = "bulk_modulus";
    const PylithReal pressureScale = _normalizer->pressureScale();

    pylith::topology::Field::Description description;
    description.label = fieldName;
    description.vectorFieldType = pylith::topology::Field::SCALAR;
    description.numComponents = 1;
    description.componentNames.resize(1);
    description.componentNames[0] = fieldName;
    description.scale = pressureScale;
    description.validator = NULL;

    _field->subfieldAdd(description, subfieldDiscretization(fieldName));
    _subfieldQueryFn(fieldName, pylith::materials::Query::dbQueryBulkModulus);

    PYLITH_METHOD_END;
} // bulkModulus


// ----------------------------------------------------------------------
// Add gravity subfield to auxiliary fields.
void
pylith::materials::AuxiliaryFactory::gravityField(spatialdata::spatialdb::GravityField* gf)
{ // gravityField
    PYLITH_METHOD_BEGIN;
    PYLITH_JOURNAL_DEBUG("gravityField(void)");

    const char* fieldName = "gravitational_acceleration";
    const char* componentNames[3] = { "gravitational_acceleration_x", "gravitational_acceleration_y", "gravitational_acceleration_z" };

    const PylithReal lengthScale = _normalizer->lengthScale();
    const PylithReal timeScale = _normalizer->timeScale();
    const PylithReal accelerationScale = lengthScale / (timeScale * timeScale);

    pylith::topology::Field::Description description;
    description.label = fieldName;
    description.vectorFieldType = pylith::topology::Field::VECTOR;
    description.numComponents = _spaceDim;
    description.componentNames.resize(_spaceDim);
    for (int i = 0; i < _spaceDim; ++i) {
        description.componentNames[i] = componentNames[i];
    } // for
    description.scale = accelerationScale;
    description.validator = NULL;

    _field->subfieldAdd(description, subfieldDiscretization(fieldName));
    _subfieldQueryFn(fieldName, pylith::materials::Query::dbQueryGravityField, gf);

    PYLITH_METHOD_END;
} // gravityField

// ----------------------------------------------------------------------
// Add body force subfield to auxiliary fields.
void
pylith::materials::AuxiliaryFactory::bodyForce(void)
{ // bodyForce
    PYLITH_METHOD_BEGIN;
    PYLITH_JOURNAL_DEBUG("bodyForce(void)");

    const char* fieldName = "body_force";
    const char* componentNames[3] = { "body_force_x", "body_force_y", "body_force_z" };

    const PylithReal forceScale = _normalizer->pressureScale() / _normalizer->lengthScale();

    pylith::topology::Field::Description description;
    description.label = fieldName;
    description.vectorFieldType = pylith::topology::Field::VECTOR;
    description.numComponents = _spaceDim;
    description.componentNames.resize(_spaceDim);
    for (int i = 0; i < _spaceDim; ++i) {
        description.componentNames[i] = componentNames[i];
    } // for
    description.scale = forceScale;
    description.validator = NULL;

    _field->subfieldAdd(description, subfieldDiscretization(fieldName));
    _subfieldQueryFn(fieldName, pylith::topology::FieldQuery::dbQueryGeneric);

    PYLITH_METHOD_END;
} // if


// ----------------------------------------------------------------------
// Add reference stress subfield to auxiliary fields.
void
pylith::materials::AuxiliaryFactory::referenceStress(void)
{ // referenceStress
    PYLITH_METHOD_BEGIN;
    PYLITH_JOURNAL_DEBUG("referenceStress(void)");

    const char* fieldName = "reference_stress";
    const char* componentNames[6] = { "reference_stress_xx", "reference_stress_yy", "reference_stress_zz", "reference_stress_xy", "reference_stress_yz", "reference_stress_xz" };
    const int stressSize = (3 == _spaceDim) ? 6 : (2 == _spaceDim) ? 4 : 1;
    const PylithReal pressureScale = _normalizer->pressureScale();

    pylith::topology::Field::Description description;
    description.label = fieldName;
    description.vectorFieldType = pylith::topology::Field::OTHER;
    description.numComponents = stressSize;
    description.componentNames.resize(stressSize);
    for (int i = 0; i < stressSize; ++i) {
        description.componentNames[i] = componentNames[i];
    } // for
    description.scale = pressureScale;
    description.validator = NULL;

    _field->subfieldAdd(description, subfieldDiscretization(fieldName));
    _subfieldQueryFn(fieldName, pylith::topology::FieldQuery::dbQueryGeneric);

    PYLITH_METHOD_END;
} // referenceStress


// ----------------------------------------------------------------------
// Add reference strain subfield to auxiliary fields.
void
pylith::materials::AuxiliaryFactory::referenceStrain(void)
{ // referenceStrain
    PYLITH_METHOD_BEGIN;
    PYLITH_JOURNAL_DEBUG("refrenceStrain(void)");

    const char* fieldName = "reference_strain";
    const char* componentNames[6] = { "reference_strain_xx", "reference_strain_yy", "reference_strain_zz", "reference_strain_xy", "reference_strain_yz", "reference_strain_xz" };
    const int strainSize = (3 == _spaceDim) ? 6 : (2 == _spaceDim) ? 4 : 1;

    pylith::topology::Field::Description description;
    description.label = fieldName;
    description.vectorFieldType = pylith::topology::Field::OTHER;
    description.numComponents = strainSize;
    description.componentNames.resize(strainSize);
    for (int i = 0; i < strainSize; ++i) {
        description.componentNames[i] = componentNames[i];
    } // for
    description.scale = 1.0;
    description.validator = NULL;

    _field->subfieldAdd(description, subfieldDiscretization(fieldName));
    _subfieldQueryFn(fieldName, pylith::topology::FieldQuery::dbQueryGeneric);

    PYLITH_METHOD_END;
} // referenceStrain


// ----------------------------------------------------------------------
// Add Maxwell time subfield to auxiliary fields.
void
pylith::materials::AuxiliaryFactory::maxwellTime(const char* identifier)
{ // maxwellTime
    PYLITH_METHOD_BEGIN;
    PYLITH_JOURNAL_DEBUG("maxwellTime(identifier)");

    const char* fieldName = "maxwell_time";
    const std::string& fieldNameFull = (identifier) ? std::string(fieldName) + std::string("_") + std::string(identifier) : fieldName;
    const PylithReal timeScale = _normalizer->timeScale();

    pylith::topology::Field::Description description;
    description.label = fieldNameFull;
    description.vectorFieldType = pylith::topology::Field::SCALAR;
    description.numComponents = 1;
    description.componentNames.resize(1);
    description.componentNames[0] = fieldNameFull;
    description.scale = timeScale;
    description.validator = pylith::topology::FieldQuery::validatorPositive;

    _field->subfieldAdd(description, subfieldDiscretization(fieldName));
    _subfieldQueryFn(fieldName, pylith::materials::Query::dbQueryMaxwellTime);

    PYLITH_METHOD_END;
} // maxwellTime

// ----------------------------------------------------------------------
// Add total strain subfield to auxiliary fields.
void
pylith::materials::AuxiliaryFactory::totalStrain(void)
{ // totalStrain
    PYLITH_METHOD_BEGIN;
    PYLITH_JOURNAL_DEBUG("totalStrain(void)");

    const char* fieldName = "total_strain";
    const char* componentNames[6] = { "total_strain_xx", "total_strain_yy", "total_strain_zz", "total_strain_xy", "total_strain_yz", "total_strain_xz" };
    const int strainSize = (3 == _spaceDim) ? 6 : (2 == _spaceDim) ? 4 : 1;

    pylith::topology::Field::Description description;
    description.label = fieldName;
    description.vectorFieldType = pylith::topology::Field::OTHER;
    description.numComponents = strainSize;
    description.componentNames.resize(strainSize);
    for (int i = 0; i < strainSize; ++i) {
        description.componentNames[i] = componentNames[i];
    } // for
    description.scale = 1.0;
    description.validator = NULL;

    _field->subfieldAdd(description, subfieldDiscretization(fieldName));
    _subfieldQueryFn(fieldName, pylith::topology::FieldQuery::dbQueryGeneric);

    PYLITH_METHOD_END;
} // totalStrain

// ----------------------------------------------------------------------
// Add viscous strain subfield to auxiliary fields.
void
pylith::materials::AuxiliaryFactory::viscousStrain(const char* identifier)
{ // viscousStrain
    PYLITH_METHOD_BEGIN;
    PYLITH_JOURNAL_DEBUG("viscousStrain(identifier)");

    const char* fieldName = "viscous_strain";
    const char* componentSuffixes[6] = { "_xx", "_yy", "_zz", "_xy", "_yz", "_xz" };
    const std::string& fieldNameFull = (identifier) ? std::string(fieldName) + std::string("_") + std::string(identifier) : fieldName;

    const int strainSize = (3 == _spaceDim) ? 6 : (2 == _spaceDim) ? 4 : 1;

    pylith::topology::Field::Description description;
    description.label = fieldName;
    description.vectorFieldType = pylith::topology::Field::OTHER;
    description.numComponents = strainSize;
    description.componentNames.resize(strainSize);
    for (int i = 0; i < strainSize; ++i) {
      description.componentNames[i] = std::string(fieldNameFull) + std::string(componentSuffixes[i]);
    } // for
    description.scale = 1.0;
    description.validator = NULL;

    _field->subfieldAdd(description, subfieldDiscretization(fieldName));
    _subfieldQueryFn(fieldName, pylith::topology::FieldQuery::dbQueryGeneric);

    PYLITH_METHOD_END;
} // viscousStrain

// ----------------------------------------------------------------------
// Add shear modulus ratio subfield to auxiliary fields.
void
pylith::materials::AuxiliaryFactory::shearModulusRatio(const char* identifier)
{ // shearModulusRatio
    PYLITH_METHOD_BEGIN;
    PYLITH_JOURNAL_DEBUG("shearModulusRatio(identifier)");

    const char* fieldName = "shear_modulus_ratio";
    const std::string& fieldNameFull = (identifier) ? std::string(fieldName) + std::string("_") + std::string(identifier) : fieldName;

    pylith::topology::Field::Description description;
    description.label = fieldNameFull;
    description.vectorFieldType = pylith::topology::Field::SCALAR;
    description.numComponents = 1;
    description.componentNames.resize(1);
    description.componentNames[0] = fieldNameFull;
    description.scale = 1.0;
    description.validator = NULL;

    _field->subfieldAdd(description, subfieldDiscretization(fieldName));
    _subfieldQueryFn(fieldName, pylith::topology::FieldQuery::dbQueryGeneric);

    PYLITH_METHOD_END;
} // shearModulusRatio

// End of file
