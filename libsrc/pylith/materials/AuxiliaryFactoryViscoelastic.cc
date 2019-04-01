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

#include "AuxiliaryFactoryViscoelastic.hh" // implementation of object methods

#include "Material.hh" // USES Material
#include "Query.hh" // USES Query

#include "pylith/topology/Field.hh" // USES Field
#include "pylith/topology/FieldQuery.hh" // HOLDSA FieldQuery

#include "spatialdata/units/Nondimensional.hh" // USES Nondimensional
#include "spatialdata/spatialdb/GravityField.hh" // USES GravityField

#include "pylith/utils/journals.hh" // USES PYLITH_JOURNAL*

#include <cassert>

// ----------------------------------------------------------------------
// Default constructor.
pylith::materials::AuxiliaryFactoryViscoelastic::AuxiliaryFactoryViscoelastic(void) {
    GenericComponent::setName("auxiliaryfactoryviscoelastic");
} // constructor


// ----------------------------------------------------------------------
// Destructor.
pylith::materials::AuxiliaryFactoryViscoelastic::~AuxiliaryFactoryViscoelastic(void) {}


// ---------------------------------------------------------------------------------------------------------------------
// Add Maxwell time subfield to auxiliary fields.
void
pylith::materials::AuxiliaryFactoryViscoelastic::addMaxwellTime(void) {
    PYLITH_METHOD_BEGIN;
    PYLITH_JOURNAL_DEBUG("addMaxwellTime(void)");

    const char* fieldName = "maxwell_time";
    const PylithReal timeScale = _normalizer->timeScale();

    pylith::topology::Field::Description description;
    description.label = fieldName;
    description.alias = fieldName;
    description.vectorFieldType = pylith::topology::Field::SCALAR;
    description.numComponents = 1;
    description.componentNames.resize(1);
    description.componentNames[0] = fieldName;
    description.scale = timeScale;
    description.validator = pylith::topology::FieldQuery::validatorPositive;

    _field->subfieldAdd(description, getSubfieldDiscretization(fieldName));
    _setSubfieldQueryFn(fieldName, pylith::materials::Query::dbQueryMaxwellTime);

    PYLITH_METHOD_END;
} // addMaxwellTime


// ---------------------------------------------------------------------------------------------------------------------
// Add Maxwell time subfield for Generalized Maxwell model to auxiliary fields.
void
pylith::materials::AuxiliaryFactoryViscoelastic::addMaxwellTimeGeneralizedMaxwell(void) {
    PYLITH_METHOD_BEGIN;
    PYLITH_JOURNAL_DEBUG("addMaxwellTimeGeneralizedMaxwell(void)");

    const char* fieldName = "maxwell_time";
    const char* componentNames[3] = { "maxwell_time_1", "maxwell_time_2", "maxwell_time_3" };
    const PylithReal timeScale = _normalizer->timeScale();

    pylith::topology::Field::Description description;
    description.label = fieldName;
    description.alias = fieldName;
    description.vectorFieldType = pylith::topology::Field::OTHER;
    description.numComponents = 3;
    description.componentNames.resize(3);
    for (int i = 0; i < 3; ++i) {
        description.componentNames[i] = componentNames[i];
    } // for
    description.scale = timeScale;
    description.validator = pylith::topology::FieldQuery::validatorPositive;

    _field->subfieldAdd(description, getSubfieldDiscretization(fieldName));
    _setSubfieldQueryFn(fieldName, pylith::materials::Query::dbQueryMaxwellTimeGeneralizedMaxwell);

    PYLITH_METHOD_END;
} // addMaxwellTimeGeneralizedMaxwell


// ---------------------------------------------------------------------------------------------------------------------
// Add shear modulus ratio subfield for generalized Maxwell model to auxiliary fields.
void
pylith::materials::AuxiliaryFactoryViscoelastic::addShearModulusRatioGeneralizedMaxwell(void) {
    PYLITH_METHOD_BEGIN;
    PYLITH_JOURNAL_DEBUG("addShearModulusRatioGeneralizedMaxwell(void)");

    const char* fieldName = "shear_modulus_ratio";

    pylith::topology::Field::Description description;
    description.label = fieldName;
    description.alias = fieldName;
    description.vectorFieldType = pylith::topology::Field::OTHER;
    description.numComponents = 3;
    description.componentNames.resize(3);
    const char* componentNames[3] = { "shear_modulus_ratio_1", "shear_modulus_ratio_2",
                                      "shear_modulus_ratio_3" };
    for (int i = 0; i < 3; ++i) {
        description.componentNames[i] = componentNames[i];
    } // for
    description.scale = 1.0;
    description.validator = NULL;

    _field->subfieldAdd(description, getSubfieldDiscretization(fieldName));
    _setSubfieldQueryFn(fieldName, pylith::materials::Query::dbQueryShearModulusRatioGeneralizedMaxwell);

    PYLITH_METHOD_END;
} // addShearModulusRatioGeneralizedMaxwell


// ---------------------------------------------------------------------------------------------------------------------
// Add total strain subfield to auxiliary fields.
void
pylith::materials::AuxiliaryFactoryViscoelastic::addTotalStrain(void) {
    PYLITH_METHOD_BEGIN;
    PYLITH_JOURNAL_DEBUG("addTotalStrain(void)");

    const char* fieldName = "total_strain";
    const char* componentNames[6] = { "total_strain_xx", "total_strain_yy", "total_strain_zz", "total_strain_xy", "total_strain_yz", "total_strain_xz" };
    const int strainSize = (3 == _spaceDim) ? 6 : (2 == _spaceDim) ? 4 : 1;

    pylith::topology::Field::Description description;
    description.label = fieldName;
    description.alias = fieldName;
    description.vectorFieldType = pylith::topology::Field::OTHER;
    description.numComponents = strainSize;
    description.componentNames.resize(strainSize);
    description.hasHistory = true;
    description.historySize = 1;
    for (int i = 0; i < strainSize; ++i) {
        description.componentNames[i] = componentNames[i];
    } // for
    description.scale = 1.0;
    description.validator = NULL;

    _field->subfieldAdd(description, getSubfieldDiscretization(fieldName));
    _setSubfieldQueryFn(fieldName, pylith::topology::FieldQuery::dbQueryGeneric);

    PYLITH_METHOD_END;
} // addTotalStrain


// ---------------------------------------------------------------------------------------------------------------------
// Add viscous strain subfield to auxiliary fields.
void
pylith::materials::AuxiliaryFactoryViscoelastic::addViscousStrain(void) {
    PYLITH_METHOD_BEGIN;
    PYLITH_JOURNAL_DEBUG("addViscousStrain(void)");

    const char* fieldName = "viscous_strain";
    const char* componentSuffixes[6] = { "_xx", "_yy", "_zz", "_xy", "_yz", "_xz" };
    const int strainSize = (3 == _spaceDim) ? 6 : (2 == _spaceDim) ? 4 : 1;

    pylith::topology::Field::Description description;
    description.label = fieldName;
    description.alias = fieldName;
    description.vectorFieldType = pylith::topology::Field::OTHER;
    description.numComponents = strainSize;
    description.componentNames.resize(strainSize);
    description.hasHistory = true;
    description.historySize = 1;
    for (int i = 0; i < strainSize; ++i) {
        description.componentNames[i] = std::string(fieldName) + std::string(componentSuffixes[i]);
    } // for
    description.scale = 1.0;
    description.validator = NULL;

    _field->subfieldAdd(description, getSubfieldDiscretization(fieldName));
    _setSubfieldQueryFn(fieldName, pylith::topology::FieldQuery::dbQueryGeneric);

    PYLITH_METHOD_END;
} // addViscousStrain


// ---------------------------------------------------------------------------------------------------------------------
// Add Generalized Maxwell viscous strain subfield to auxiliary fields.
void
pylith::materials::AuxiliaryFactoryViscoelastic::addViscousStrainGeneralizedMaxwell(void) {
    PYLITH_METHOD_BEGIN;
    PYLITH_JOURNAL_DEBUG("addViscousStrainGeneralizedMaxwell(void)");

    const char* fieldName = "viscous_strain";
    const char* componentElementNumbers[3] = { "_1", "_2", "_3" };
    const char* componentSuffixes[6] = { "_xx", "_yy", "_zz", "_xy", "_yz", "_xz" };
    const int strainSize = (3 == _spaceDim) ? 6 : (2 == _spaceDim) ? 4 : 1;

    pylith::topology::Field::Description description;
    description.label = fieldName;
    description.alias = fieldName;
    description.vectorFieldType = pylith::topology::Field::OTHER;
    description.numComponents = 3 * strainSize;
    description.componentNames.resize(3 * strainSize);
    description.hasHistory = true;
    description.historySize = 1;
    for (int j = 0, iname = 0; j < 3; ++j) {
        for (int i = 0; i < strainSize; ++i, ++iname) {
            description.componentNames[iname] = std::string(fieldName) + std::string(componentElementNumbers[j]) + std::string(componentSuffixes[i]);
        } // for i
    } // for j
    description.scale = 1.0;
    description.validator = NULL;

    _field->subfieldAdd(description, getSubfieldDiscretization(fieldName));
    _setSubfieldQueryFn(fieldName, pylith::topology::FieldQuery::dbQueryGeneric);

    PYLITH_METHOD_END;
} // addViscousStrainGeneralizedMaxwell


// End of file
