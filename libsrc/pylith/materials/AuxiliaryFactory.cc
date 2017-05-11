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

#include "MaterialNew.hh" // USES MaterialNew
#include "Query.hh" // USES Query

#include "pylith/topology/Field.hh" // USES Field
#include "pylith/topology/FieldQuery.hh" // HOLDSA FieldQuery

#include "spatialdata/units/Nondimensional.hh" // USES Nondimensional
#include "spatialdata/spatialdb/GravityField.hh" // USES GravityField

#include "pylith/utils/journals.hh" // USES PYLITH_JOURNAL*

#include <cassert>

// ----------------------------------------------------------------------
const char* pylith::materials::AuxiliaryFactory::_genericComponent = "auxiliaryfactory";

// ----------------------------------------------------------------------
// Default constructor.
pylith::materials::AuxiliaryFactory::AuxiliaryFactory(const MaterialNew& material,
                                                      const spatialdata::units::Nondimensional& normalizer,
                                                      const int spaceDim) :
    _material(material),
    _normalizer(normalizer),
    _spaceDim(spaceDim)
{ // constructor
    assert(1 <= _spaceDim && _spaceDim <= 3);
    GenericComponent::name(_genericComponent);
} // constructor

// ----------------------------------------------------------------------
// Destructor.
pylith::materials::AuxiliaryFactory::~AuxiliaryFactory(void)
{ // destructor
} // destructor

// ----------------------------------------------------------------------
// Add density field to auxiliary fields.
void
pylith::materials::AuxiliaryFactory::density(void) const
{ // density
    PYLITH_METHOD_BEGIN;
    PYLITH_JOURNAL_DEBUG("density(void)");

    const char* fieldName = "density";
    const PylithReal densityScale = _normalizer.densityScale();

    pylith::topology::Field::Description description;
    description.label = fieldName;
    description.vectorFieldType = pylith::topology::Field::SCALAR;
    description.numComponents = 1;
    description.componentNames.resize(1);
    description.componentNames[0] = fieldName;
    description.scale = densityScale;
    description.validator = pylith::topology::FieldQuery::validatorPositive;

    assert(_material._auxFields);
    assert(_material._auxFieldsQuery);
    _material._auxFields->subfieldAdd(description, _material.auxFieldDiscretization(fieldName));
    _material._auxFieldsQuery->queryFn(fieldName, pylith::topology::FieldQuery::dbQueryGeneric);

    PYLITH_METHOD_END;
} // density


// ----------------------------------------------------------------------
// Add shear modulus field to auxiliary fields.
void
pylith::materials::AuxiliaryFactory::shearModulus(void) const
{ // shearModulus
    PYLITH_METHOD_BEGIN;
    PYLITH_JOURNAL_DEBUG("shearModulus(void)");

    const char* fieldName = "shear_modulus";
    const PylithReal pressureScale = _normalizer.pressureScale();

    pylith::topology::Field::Description description;
    description.label = fieldName;
    description.vectorFieldType = pylith::topology::Field::SCALAR;
    description.numComponents = 1;
    description.componentNames.resize(1);
    description.componentNames[0] = fieldName;
    description.scale = pressureScale;
    description.validator = NULL;

    assert(_material._auxFields);
    assert(_material._auxFieldsQuery);
    _material._auxFields->subfieldAdd(description, _material.auxFieldDiscretization(fieldName));
    _material._auxFieldsQuery->queryFn(fieldName, pylith::materials::Query::dbQueryShearModulus);

    PYLITH_METHOD_END;
} // shearModulus


// ----------------------------------------------------------------------
// Add bulk modulus field to auxiliary fields.
void
pylith::materials::AuxiliaryFactory::bulkModulus(void) const
{ // bulkModulus
    PYLITH_METHOD_BEGIN;
    PYLITH_JOURNAL_DEBUG("bulkModulus(void)");

    const char* fieldName = "bulk_modulus";
    const PylithReal pressureScale = _normalizer.pressureScale();

    pylith::topology::Field::Description description;
    description.label = fieldName;
    description.vectorFieldType = pylith::topology::Field::SCALAR;
    description.numComponents = 1;
    description.componentNames.resize(1);
    description.componentNames[0] = fieldName;
    description.scale = pressureScale;
    description.validator = NULL;

    assert(_material._auxFields);
    assert(_material._auxFieldsQuery);
    _material._auxFields->subfieldAdd(description, _material.auxFieldDiscretization(fieldName));
    _material._auxFieldsQuery->queryFn(fieldName, pylith::materials::Query::dbQueryBulkModulus);

    PYLITH_METHOD_END;
} // bulkModulus


// ----------------------------------------------------------------------
// Add gravity field to auxiliary fields.
void
pylith::materials::AuxiliaryFactory::gravityField(spatialdata::spatialdb::GravityField* gf) const
{ // gravityField
    PYLITH_METHOD_BEGIN;
    PYLITH_JOURNAL_DEBUG("gravityField(void)");

    const char* fieldName = "gravity_field";
    const char* componentNames[3] = { "gravity_field_x", "gravity_field_y", "gravity_field_z" };

    const PylithReal lengthScale = _normalizer.lengthScale();
    const PylithReal timeScale = _normalizer.timeScale();
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

    assert(_material._auxFields);
    assert(_material._auxFieldsQuery);
    _material._auxFields->subfieldAdd(description, _material.auxFieldDiscretization(fieldName));
    _material._auxFieldsQuery->queryFn(fieldName, pylith::materials::Query::dbQueryGravityField, gf);

    PYLITH_METHOD_END;
} // gravityField

// ----------------------------------------------------------------------
// Add body force field to auxiliary fields.
void
pylith::materials::AuxiliaryFactory::bodyForce(void) const
{ // bodyForce
    PYLITH_METHOD_BEGIN;
    PYLITH_JOURNAL_DEBUG("bodyForce(void)");

    const char* fieldName = "body_force";
    const char* componentNames[3] = { "body_force_x", "body_force_y", "body_force_z" };

    const PylithReal densityScale = _normalizer.densityScale();
    const PylithReal lengthScale = _normalizer.lengthScale();
    const PylithReal timeScale = _normalizer.timeScale();
    const PylithReal accelerationScale = lengthScale / (timeScale * timeScale);
    const PylithReal forceScale = densityScale * accelerationScale;

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

    assert(_material._auxFields);
    assert(_material._auxFieldsQuery);
    _material._auxFields->subfieldAdd(description, _material.auxFieldDiscretization(fieldName));
    _material._auxFieldsQuery->queryFn(fieldName, pylith::topology::FieldQuery::dbQueryGeneric);

    PYLITH_METHOD_END;
} // if


// ----------------------------------------------------------------------
// Add reference stress to auxiliary fields.
void
pylith::materials::AuxiliaryFactory::referenceStress(void) const
{ // referenceStress
    PYLITH_METHOD_BEGIN;
    PYLITH_JOURNAL_DEBUG("referenceStress(void)");

    const char* fieldName = "reference_stress";
    const char* componentNames[6] = { "stress_xx", "stress_yy", "stress_zz", "stress_xy", "stress_yz", "stress_xz" };
    const int stressSize = (3 == _spaceDim) ? 6 : (2 == _spaceDim) ? 4 : 1;
    const PylithReal pressureScale = _normalizer.pressureScale();

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

    assert(_material._auxFields);
    assert(_material._auxFieldsQuery);
    _material._auxFields->subfieldAdd(description, _material.auxFieldDiscretization(fieldName));
    _material._auxFieldsQuery->queryFn(fieldName, pylith::topology::FieldQuery::dbQueryGeneric);

    PYLITH_METHOD_END;
} // referenceStress


// ----------------------------------------------------------------------
// Add reference strain to auxiliary fields.
void
pylith::materials::AuxiliaryFactory::referenceStrain(void) const
{ // referenceStrain
    PYLITH_METHOD_BEGIN;
    PYLITH_JOURNAL_DEBUG("refrenceStrain(void)");

    const char* fieldName = "reference_strain";
    const char* componentNames[6] = { "strain_xx", "strain_yy", "strain_zz", "strain_xy", "strain_yz", "strain_xz" };
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

    assert(_material._auxFields);
    assert(_material._auxFieldsQuery);
    _material._auxFields->subfieldAdd(description, _material.auxFieldDiscretization(fieldName));
    _material._auxFieldsQuery->queryFn(fieldName, pylith::topology::FieldQuery::dbQueryGeneric);

    PYLITH_METHOD_END;
} // referenceStrain


// ----------------------------------------------------------------------
// Add Maxwell time field to auxiliary fields.
void
pylith::materials::AuxiliaryFactory::maxwellTime(void) const
{ // maxwellTime
    PYLITH_METHOD_BEGIN;
    PYLITH_JOURNAL_DEBUG("maxwellTime(void)");

    const char* fieldName = "maxwell_time";
    const PylithReal timeScale = _normalizer.timeScale();

    pylith::topology::Field::Description description;
    description.label = fieldName;
    description.vectorFieldType = pylith::topology::Field::SCALAR;
    description.numComponents = 1;
    description.componentNames.resize(1);
    description.componentNames[0] = fieldName;
    description.scale = timeScale;
    description.validator = pylith::topology::FieldQuery::validatorPositive;

    assert(_material._auxFields);
    assert(_material._auxFieldsQuery);
    _material._auxFields->subfieldAdd(description, _material.auxFieldDiscretization(fieldName));
    _material._auxFieldsQuery->queryFn(fieldName, pylith::topology::FieldQuery::dbQueryGeneric);

    PYLITH_METHOD_END;
} // maxwellTime

// ----------------------------------------------------------------------
// Add total strain field to auxiliary fields.
void
pylith::materials::AuxiliaryFactory::totalStrain(void) const
{ // totalStrain
    PYLITH_METHOD_BEGIN;
    PYLITH_JOURNAL_DEBUG("totalStrain(void)");

    const char* fieldName = "total_strain";
    const char* componentNames[6] = { "strain_xx", "strain_yy", "strain_zz", "strain_xy", "strain_yz", "strain_xz" };
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

    assert(_material._auxFields);
    assert(_material._auxFieldsQuery);
    _material._auxFields->subfieldAdd(description, _material.auxFieldDiscretization(fieldName));
    _material._auxFieldsQuery->queryFn(fieldName, pylith::topology::FieldQuery::dbQueryGeneric);

    PYLITH_METHOD_END;
} // totalStrain

// ----------------------------------------------------------------------
// Add viscous strain field to auxiliary fields.
void
pylith::materials::AuxiliaryFactory::viscousStrain(void) const
{ // viscousStrain
    PYLITH_METHOD_BEGIN;
    PYLITH_JOURNAL_DEBUG("viscousStrain(void)");

    const char* fieldName = "viscous_strain";
    const char* componentNames[6] = { "strain_xx", "strain_yy", "strain_zz", "strain_xy", "strain_yz", "strain_xz" };
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

    assert(_material._auxFields);
    assert(_material._auxFieldsQuery);
    _material._auxFields->subfieldAdd(description, _material.auxFieldDiscretization(fieldName));
    _material._auxFieldsQuery->queryFn(fieldName, pylith::topology::FieldQuery::dbQueryGeneric);

    PYLITH_METHOD_END;
} // viscousStrain

// End of file
