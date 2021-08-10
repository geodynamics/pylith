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

#include "AuxiliaryFactoryPoroelasticity.hh" // implementation of object methods

#include "Query.hh" // USES Query

#include "pylith/topology/Field.hh" // USES Field
#include "pylith/topology/FieldQuery.hh" // HOLDSA FieldQuery

#include "spatialdata/units/Nondimensional.hh" // USES Nondimensional
#include "spatialdata/spatialdb/GravityField.hh" // USES GravityField

#include "pylith/utils/error.hh" // USES PYLITH_METHOD*
#include "pylith/utils/journals.hh" // USES PYLITH_JOURNAL*

#include <cassert>

// ---------------------------------------------------------------------------------------------------------------------
// Default constructor.
pylith::materials::AuxiliaryFactoryPoroelasticity::AuxiliaryFactoryPoroelasticity(void) {
    GenericComponent::setName("auxiliaryfactoryporoelasticity");
} // constructor


// ---------------------------------------------------------------------------------------------------------------------
// Destructor.
pylith::materials::AuxiliaryFactoryPoroelasticity::~AuxiliaryFactoryPoroelasticity(void) {}


// ---------------------------------------------------------------------------------------------------------------------
// Add body force subfield to auxiliary fields.
void
pylith::materials::AuxiliaryFactoryPoroelasticity::addBodyForce(void) {
    PYLITH_METHOD_BEGIN;
    PYLITH_JOURNAL_DEBUG("addBodyForce(void)");

    const char *subfieldName = "body_force";
    const char *componentNames[3] = {
        "body_force_x",
        "body_force_y",
        "body_force_z"
    };

    const PylithReal forceScale = _normalizer->getPressureScale() / _normalizer->getLengthScale();

    pylith::topology::Field::Description description;
    description.label = subfieldName;
    description.alias = subfieldName;
    description.vectorFieldType = pylith::topology::Field::VECTOR;
    description.numComponents = _spaceDim;
    description.componentNames.resize(_spaceDim);
    for (int i = 0; i < _spaceDim; ++i) {
        description.componentNames[i] = componentNames[i];
    } // for
    description.scale = forceScale;
    description.validator = NULL;

    _field->subfieldAdd(description, getSubfieldDiscretization(subfieldName));
    this->setSubfieldQuery(subfieldName);

    PYLITH_METHOD_END;
} // addBodyForce


// ---------------------------------------------------------------------------------------------------------------------
// Add gravity subfield to auxiliary fields.
void
pylith::materials::AuxiliaryFactoryPoroelasticity::addGravityField(spatialdata::spatialdb::GravityField *gf) {
    PYLITH_METHOD_BEGIN;
    PYLITH_JOURNAL_DEBUG("addGravityField(void)");

    const char *subfieldName = "gravitational_acceleration";
    const char *componentNames[3] = {
        "gravitational_acceleration_x",
        "gravitational_acceleration_y",
        "gravitational_acceleration_z"
    };
    const PylithReal lengthScale = _normalizer->getLengthScale();
    const PylithReal timeScale = _normalizer->getTimeScale();
    const PylithReal accelerationScale = lengthScale / (timeScale * timeScale);

    pylith::topology::Field::Description description;
    description.label = subfieldName;
    description.alias = subfieldName;
    description.vectorFieldType = pylith::topology::Field::VECTOR;
    description.numComponents = _spaceDim;
    description.componentNames.resize(_spaceDim);
    for (int i = 0; i < _spaceDim; ++i) {
        description.componentNames[i] = componentNames[i];
    } // for
    description.scale = accelerationScale;
    description.validator = NULL;

    _field->subfieldAdd(description, getSubfieldDiscretization(subfieldName));
    pylith::materials::Query::gravityFieldFromDB(subfieldName, this, gf, _spaceDim);

    PYLITH_METHOD_END;
} // addGravityField


// ----------------------------------------------------------------------
// Add porosity subfield to auxiliary fields.
void
pylith::materials::AuxiliaryFactoryPoroelasticity::addPorosity(void) { // porosity
    PYLITH_METHOD_BEGIN;
    PYLITH_JOURNAL_DEBUG("addPorosity(void)");

    const char *subfieldName = "porosity";
    const PylithReal noScale = 1;

    pylith::topology::Field::Description description;
    description.label = subfieldName;
    description.alias = subfieldName;
    description.vectorFieldType = pylith::topology::Field::SCALAR;
    description.numComponents = 1;
    description.componentNames.resize(1);
    description.hasHistory = true;
    description.historySize = 1;
    description.componentNames[0] = subfieldName;
    description.scale = noScale;
    description.validator = pylith::topology::FieldQuery::validatorNonnegative;

    _field->subfieldAdd(description, getSubfieldDiscretization(subfieldName));
    this->setSubfieldQuery(subfieldName);

    PYLITH_METHOD_END;
} // addPorosity


// ---------------------------------------------------------------------------------------------------------------------
// Add density subfield to auxiliary fields.
void
pylith::materials::AuxiliaryFactoryPoroelasticity::addSolidDensity(void) {
    PYLITH_METHOD_BEGIN;
    PYLITH_JOURNAL_DEBUG("addSolidDensity(void)");

    const char *subfieldName = "solid_density";
    const PylithReal densityScale = _normalizer->getDensityScale();

    pylith::topology::Field::Description description;
    description.label = subfieldName;
    description.alias = subfieldName;
    description.vectorFieldType = pylith::topology::Field::SCALAR;
    description.numComponents = 1;
    description.componentNames.resize(1);
    description.componentNames[0] = subfieldName;
    description.scale = densityScale;
    description.validator = pylith::topology::FieldQuery::validatorPositive;

    _field->subfieldAdd(description, getSubfieldDiscretization(subfieldName));
    this->setSubfieldQuery(subfieldName);

    PYLITH_METHOD_END;
} // addSolidDensity


// ----------------------------------------------------------------------
// Add fluid density subfield to auxiliary fields.
void
pylith::materials::AuxiliaryFactoryPoroelasticity::addFluidDensity(void) { // fluidDensity
    PYLITH_METHOD_BEGIN;
    PYLITH_JOURNAL_DEBUG("addFluidDensity(void)");

    const char *subfieldName = "fluid_density";
    const PylithReal densityScale = _normalizer->getDensityScale();

    pylith::topology::Field::Description description;
    description.label = subfieldName;
    description.alias = subfieldName;
    description.vectorFieldType = pylith::topology::Field::SCALAR;
    description.numComponents = 1;
    description.componentNames.resize(1);
    description.componentNames[0] = subfieldName;
    description.scale = densityScale;
    description.validator = pylith::topology::FieldQuery::validatorPositive;

    _field->subfieldAdd(description, getSubfieldDiscretization(subfieldName));
    this->setSubfieldQuery(subfieldName);

    PYLITH_METHOD_END;
} // addFluidDensity


// ----------------------------------------------------------------------
// Add fluid viscosity subfield to auxiliary fields.
void
pylith::materials::AuxiliaryFactoryPoroelasticity::addFluidViscosity(void) { // fluidViscosity
    PYLITH_METHOD_BEGIN;
    PYLITH_JOURNAL_DEBUG("addFluidViscosity(void)");

    const char *subfieldName = "fluid_viscosity";
    const PylithReal pressureScale = _normalizer->getPressureScale();
    const PylithReal timeScale = _normalizer->getTimeScale();
    const PylithReal viscosityScale = pressureScale * timeScale;

    pylith::topology::Field::Description description;
    description.label = subfieldName;
    description.alias = subfieldName;
    description.vectorFieldType = pylith::topology::Field::SCALAR;
    description.numComponents = 1;
    description.componentNames.resize(1);
    description.componentNames[0] = subfieldName;
    description.scale = viscosityScale;
    description.validator = pylith::topology::FieldQuery::validatorPositive;

    _field->subfieldAdd(description, getSubfieldDiscretization(subfieldName));
    this->setSubfieldQuery(subfieldName);

    PYLITH_METHOD_END;

} // addFluidViscosity


// ----------------------------------------------------------------------
// Add source density subfield to auxiliary fields.
void
pylith::materials::AuxiliaryFactoryPoroelasticity::addSourceDensity(void) { // sourceDensity
    PYLITH_METHOD_BEGIN;
    PYLITH_JOURNAL_DEBUG("addSourceDensity(void)");

    const char *subfieldName = "source_density";
    const PylithReal lengthScale = _normalizer->getLengthScale();
    const PylithReal timeScale = _normalizer->getTimeScale();
    const PylithReal sourceDensityScale = lengthScale / timeScale;

    pylith::topology::Field::Description description;
    description.label = subfieldName;
    description.alias = subfieldName;
    description.vectorFieldType = pylith::topology::Field::SCALAR;
    description.numComponents = 1;
    description.componentNames.resize(1);
    description.componentNames[0] = subfieldName;
    description.scale = sourceDensityScale;
    description.validator = pylith::topology::FieldQuery::validatorPositive;

    _field->subfieldAdd(description, getSubfieldDiscretization(subfieldName));
    this->setSubfieldQuery(subfieldName);

    PYLITH_METHOD_END;

} // addSourceDensity


// ----------------------------------------------------------------------
// Add constant pressure source subfield to auxiliary fields.
void
pylith::materials::AuxiliaryFactoryPoroelasticity::addConstantPressureSource(void) { // constantPressureSource
    PYLITH_METHOD_BEGIN;
    PYLITH_JOURNAL_DEBUG("addConstantPressureSource(void)");

    const char *subfieldName = "constant_pressure_source";
    const PylithReal pressureScale = _normalizer->getPressureScale();
    const PylithReal constantPressureSourceScale = pressureScale;

    pylith::topology::Field::Description description;
    description.label = subfieldName;
    description.alias = subfieldName;
    description.vectorFieldType = pylith::topology::Field::SCALAR;
    description.numComponents = 1;
    description.componentNames.resize(1);
    description.componentNames[0] = subfieldName;
    description.scale = constantPressureSourceScale;
    description.validator = pylith::topology::FieldQuery::validatorNonnegative;

    _field->subfieldAdd(description, getSubfieldDiscretization(subfieldName));
    this->setSubfieldQuery(subfieldName);

    PYLITH_METHOD_END;

} // addConstantPressureSource


// End of file
