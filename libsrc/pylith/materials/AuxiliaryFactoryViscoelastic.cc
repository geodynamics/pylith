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

#include "AuxiliaryFactoryViscoelastic.hh" // implementation of object methods

#include "Material.hh" // USES Material
#include "Query.hh" // USES Query

#include "pylith/topology/Field.hh" // USES Field
#include "pylith/topology/FieldQuery.hh" // HOLDSA FieldQuery

#include "spatialdata/units/Nondimensional.hh" // USES Nondimensional
#include "spatialdata/spatialdb/GravityField.hh" // USES GravityField

#include "pylith/utils/error.hh" // USES PYLITH_METHOD*
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

    const char* subfieldName = "maxwell_time";
    const PylithReal timeScale = _normalizer->getTimeScale();

    pylith::topology::Field::Description description;
    description.label = subfieldName;
    description.alias = subfieldName;
    description.vectorFieldType = pylith::topology::Field::SCALAR;
    description.numComponents = 1;
    description.componentNames.resize(1);
    description.componentNames[0] = subfieldName;
    description.scale = timeScale;
    description.validator = pylith::topology::FieldQuery::validatorPositive;

    _field->subfieldAdd(description, getSubfieldDiscretization(subfieldName));
    pylith::materials::Query::maxwellTimeFromVM(subfieldName, this);

    PYLITH_METHOD_END;
} // addMaxwellTime


// ---------------------------------------------------------------------------------------------------------------------
// Add Maxwell time subfield for Generalized Maxwell model to auxiliary fields.
void
pylith::materials::AuxiliaryFactoryViscoelastic::addMaxwellTimeGeneralizedMaxwell(void) {
    PYLITH_METHOD_BEGIN;
    PYLITH_JOURNAL_DEBUG("addMaxwellTimeGeneralizedMaxwell(void)");

    const char* subfieldName = "maxwell_time";
    const char* componentNames[3] = {
        "maxwell_time_1",
        "maxwell_time_2",
        "maxwell_time_3"
    };
    const PylithReal timeScale = _normalizer->getTimeScale();

    pylith::topology::Field::Description description;
    description.label = subfieldName;
    description.alias = subfieldName;
    description.vectorFieldType = pylith::topology::Field::OTHER;
    description.numComponents = 3;
    description.componentNames.resize(3);
    for (int i = 0; i < 3; ++i) {
        description.componentNames[i] = componentNames[i];
    } // for
    description.scale = timeScale;
    description.validator = pylith::topology::FieldQuery::validatorPositive;

    _field->subfieldAdd(description, getSubfieldDiscretization(subfieldName));
    pylith::materials::Query::generalizedMaxwellTimesFromVM(subfieldName, this);

    PYLITH_METHOD_END;
} // addMaxwellTimeGeneralizedMaxwell


// ---------------------------------------------------------------------------------------------------------------------
// Add shear modulus ratio subfield for generalized Maxwell model to auxiliary fields.
void
pylith::materials::AuxiliaryFactoryViscoelastic::addShearModulusRatioGeneralizedMaxwell(void) {
    PYLITH_METHOD_BEGIN;
    PYLITH_JOURNAL_DEBUG("addShearModulusRatioGeneralizedMaxwell(void)");

    const char* subfieldName = "shear_modulus_ratio";

    pylith::topology::Field::Description description;
    description.label = subfieldName;
    description.alias = subfieldName;
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

    _field->subfieldAdd(description, getSubfieldDiscretization(subfieldName));
    pylith::materials::Query::generalizedMaxwellShearModulusRatiosFromVM(subfieldName, this);

    PYLITH_METHOD_END;
} // addShearModulusRatioGeneralizedMaxwell


// ---------------------------------------------------------------------------------------------------------------------
// Add power-law reference strain rate subfield to auxiliary fields.
void
pylith::materials::AuxiliaryFactoryViscoelastic::addPowerLawReferenceStrainRate(void) {
    PYLITH_METHOD_BEGIN;
    PYLITH_JOURNAL_DEBUG("addPowerLawReferenceStrainRate(void)");

    const char* subfieldName = "power_law_reference_strain_rate";
    const PylithReal strainRateScale = 1.0/_normalizer->getTimeScale();

    pylith::topology::Field::Description description;
    description.label = subfieldName;
    description.alias = subfieldName;
    description.vectorFieldType = pylith::topology::Field::SCALAR;
    description.numComponents = 1;
    description.componentNames.resize(1);
    description.componentNames[0] = subfieldName;
    description.scale = strainRateScale;
    description.validator = pylith::topology::FieldQuery::validatorPositive;

    _field->subfieldAdd(description, getSubfieldDiscretization(subfieldName));
    this->setSubfieldQuery(subfieldName);

    PYLITH_METHOD_END;
} // addPowerLawReferenceStrainRate


// ---------------------------------------------------------------------------------------------------------------------
// Add power-law reference stress subfield to auxiliary fields.
void
pylith::materials::AuxiliaryFactoryViscoelastic::addPowerLawReferenceStress(void) {
    PYLITH_METHOD_BEGIN;
    PYLITH_JOURNAL_DEBUG("addPowerLawReferenceStress(void)");

    const char* subfieldName = "power_law_reference_stress";
    const PylithReal pressureScale = _normalizer->getPressureScale();

    pylith::topology::Field::Description description;
    description.label = subfieldName;
    description.alias = subfieldName;
    description.vectorFieldType = pylith::topology::Field::SCALAR;
    description.numComponents = 1;
    description.componentNames.resize(1);
    description.componentNames[0] = subfieldName;
    description.scale = pressureScale;
    description.validator = pylith::topology::FieldQuery::validatorPositive;

    _field->subfieldAdd(description, getSubfieldDiscretization(subfieldName));
    this->setSubfieldQuery(subfieldName);

    PYLITH_METHOD_END;
} // addPowerLawReferenceStress


// ---------------------------------------------------------------------------------------------------------------------
// Add power-law exponenet subfield to auxiliary fields.
void
pylith::materials::AuxiliaryFactoryViscoelastic::addPowerLawExponent(void) {
    PYLITH_METHOD_BEGIN;
    PYLITH_JOURNAL_DEBUG("addPowerLawExponent(void)");

    const char* subfieldName = "power_law_exponent";

    pylith::topology::Field::Description description;
    description.label = subfieldName;
    description.alias = subfieldName;
    description.vectorFieldType = pylith::topology::Field::SCALAR;
    description.numComponents = 1;
    description.componentNames.resize(1);
    description.componentNames[0] = subfieldName;
    description.scale = 1.0;
    description.validator = pylith::topology::FieldQuery::validatorPositive;

    _field->subfieldAdd(description, getSubfieldDiscretization(subfieldName));
    this->setSubfieldQuery(subfieldName);

    PYLITH_METHOD_END;
} // addPowerLawExponent


// ---------------------------------------------------------------------------------------------------------------------
// Add total strain subfield to auxiliary fields.
void
pylith::materials::AuxiliaryFactoryViscoelastic::addTotalStrain(void) {
    PYLITH_METHOD_BEGIN;
    PYLITH_JOURNAL_DEBUG("addTotalStrain(void)");

    const char* subfieldName = "total_strain";
    const char* componentNames[6] = {
        "total_strain_xx",
        "total_strain_yy",
        "total_strain_zz",
        "total_strain_xy",
        "total_strain_yz",
        "total_strain_xz"
    };
    const int strainSize = (3 == _spaceDim) ? 6 : (2 == _spaceDim) ? 4 : 1;

    pylith::topology::Field::Description description;
    description.label = subfieldName;
    description.alias = subfieldName;
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

    _field->subfieldAdd(description, getSubfieldDiscretization(subfieldName));
    this->setSubfieldQuery(subfieldName);

    PYLITH_METHOD_END;
} // addTotalStrain


// ---------------------------------------------------------------------------------------------------------------------
// Add stress subfield to auxiliary fields.
void
pylith::materials::AuxiliaryFactoryViscoelastic::addStress(void) {
    PYLITH_METHOD_BEGIN;
    PYLITH_JOURNAL_DEBUG("addStress(void)");

    const char* subfieldName = "stress";
    const char* componentNames[6] = {
        "stress_xx",
        "stress_yy",
        "stress_zz",
        "stress_xy",
        "stress_yz",
        "stress_xz"
    };
    const int stressSize = (3 == _spaceDim) ? 6 : (2 == _spaceDim) ? 4 : 1;

    pylith::topology::Field::Description description;
    description.label = subfieldName;
    description.alias = subfieldName;
    description.vectorFieldType = pylith::topology::Field::OTHER;
    description.numComponents = stressSize;
    description.componentNames.resize(stressSize);
    description.hasHistory = true;
    description.historySize = 1;
    for (int i = 0; i < stressSize; ++i) {
        description.componentNames[i] = componentNames[i];
    } // for
    description.scale = 1.0;
    description.validator = NULL;

    _field->subfieldAdd(description, getSubfieldDiscretization(subfieldName));
    this->setSubfieldQuery(subfieldName);

    PYLITH_METHOD_END;
} // addStress


// ---------------------------------------------------------------------------------------------------------------------
// Add viscous strain subfield to auxiliary fields.
void
pylith::materials::AuxiliaryFactoryViscoelastic::addViscousStrain(void) {
    PYLITH_METHOD_BEGIN;
    PYLITH_JOURNAL_DEBUG("addViscousStrain(void)");

    const char* subfieldName = "viscous_strain";
    const char* componentNames[6] = {
        "viscous_strain_xx",
        "viscous_strain_yy",
        "viscous_strain_zz",
        "viscous_strain_xy",
        "viscous_strain_yz",
        "viscous_strain_xz"
    };
    const int strainSize = (3 == _spaceDim) ? 6 : (2 == _spaceDim) ? 4 : 1;

    pylith::topology::Field::Description description;
    description.label = subfieldName;
    description.alias = subfieldName;
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

    _field->subfieldAdd(description, getSubfieldDiscretization(subfieldName));
    this->setSubfieldQuery(subfieldName);

    PYLITH_METHOD_END;
} // addViscousStrain


// ---------------------------------------------------------------------------------------------------------------------
// Add Generalized Maxwell viscous strain subfield to auxiliary fields.
void
pylith::materials::AuxiliaryFactoryViscoelastic::addViscousStrainGeneralizedMaxwell(void) {
    PYLITH_METHOD_BEGIN;
    PYLITH_JOURNAL_DEBUG("addViscousStrainGeneralizedMaxwell(void)");

    const char* subfieldName = "viscous_strain";
    const char* componentElementNumbers[3] = { "_1", "_2", "_3" };
    const char* componentSuffixes[6] = { "_xx", "_yy", "_zz", "_xy", "_yz", "_xz" };
    const int strainSize = (3 == _spaceDim) ? 6 : (2 == _spaceDim) ? 4 : 1;

    pylith::topology::Field::Description description;
    description.label = subfieldName;
    description.alias = subfieldName;
    description.vectorFieldType = pylith::topology::Field::OTHER;
    description.numComponents = 3 * strainSize;
    description.componentNames.resize(3 * strainSize);
    description.hasHistory = true;
    description.historySize = 1;
    for (int j = 0, iname = 0; j < 3; ++j) {
        for (int i = 0; i < strainSize; ++i, ++iname) {
            description.componentNames[iname] = std::string(subfieldName) + std::string(componentElementNumbers[j]) + std::string(componentSuffixes[i]);
        } // for i
    } // for j
    description.scale = 1.0;
    description.validator = NULL;

    _field->subfieldAdd(description, getSubfieldDiscretization(subfieldName));
    this->setSubfieldQuery(subfieldName);

    PYLITH_METHOD_END;
} // addViscousStrainGeneralizedMaxwell


// End of file
