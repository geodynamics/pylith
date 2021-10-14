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

#include "AuxiliaryFactoryWellboreSource.hh" // implementation of object methods

#include "pylith/topology/Field.hh" // USES Field
#include "pylith/topology/FieldQuery.hh" // HOLDSA FieldQuery

#include "spatialdata/units/Nondimensional.hh" // USES Nondimensional

#include "pylith/utils/journals.hh" // USES PYLITH_JOURNAL*

#include <cassert>

// ---------------------------------------------------------------------------------------------------------------------
// Default constructor.
pylith::sources::AuxiliaryFactoryWellboreSource::AuxiliaryFactoryWellboreSource(void) {
    GenericComponent::setName("auxiliaryfactorywellboresource");
} // constructor


// ---------------------------------------------------------------------------------------------------------------------
// Destructor.
pylith::sources::AuxiliaryFactoryWellboreSource::~AuxiliaryFactoryWellboreSource(void) {}


// ----------------------------------------------------------------------
// Add fluid density subfield to auxiliary fields.
void
pylith::sources::AuxiliaryFactoryWellboreSource::addFluidDensity(void) { // fluidDensity
    PYLITH_METHOD_BEGIN;
    PYLITH_JOURNAL_DEBUG("addFluidDensity(void)");

    const char* subfieldName = "fluid_density";
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
pylith::sources::AuxiliaryFactoryWellboreSource::addFluidViscosity(void) { // fluidViscosity
    PYLITH_METHOD_BEGIN;
    PYLITH_JOURNAL_DEBUG("addFluidViscosity(void)");

    const char* subfieldName = "fluid_viscosity";
    const PylithReal pressureScale = _normalizer->getPressureScale();
    const PylithReal timeScale = _normalizer->getTimeScale();
    const PylithReal viscosityScale = pressureScale*timeScale;

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


// ----------------------------------------------------------------------------
// Add isotropic permeability subfield to auxiliary fields.
void
pylith::sources::AuxiliaryFactoryWellboreSource::addIsotropicPermeability(void) { // isotropicPermeablity
    PYLITH_METHOD_BEGIN;
    PYLITH_JOURNAL_DEBUG("addIsotropicPermeability(void)");

    const char* subfieldName = "isotropic_permeability";

    const PylithReal lengthScale = _normalizer->getLengthScale();
    const PylithReal permeabilityScale = lengthScale*lengthScale;

    pylith::topology::Field::Description description;
    description.label = subfieldName;
    description.alias = subfieldName;
    description.vectorFieldType = pylith::topology::Field::SCALAR;
    description.numComponents = 1;
    description.componentNames.resize(1);
    description.componentNames[0] = subfieldName;
    description.scale = permeabilityScale;
    description.validator = NULL;

    _field->subfieldAdd(description, getSubfieldDiscretization(subfieldName));
    this->setSubfieldQuery(subfieldName);

    PYLITH_METHOD_END;
} // addIsotropicPermeability


// ----------------------------------------------------------------------------
// Add wellbore radius subfield to auxiliary fields.
void
pylith::sources::AuxiliaryFactoryWellboreSource::addWellboreRadius(void) { // wellboreRadius
    PYLITH_METHOD_BEGIN;
    PYLITH_JOURNAL_DEBUG("addWellboreRadius(void)");

    const char* subfieldName = "wellbore_radius";

    const PylithReal lengthScale = _normalizer->getLengthScale();

    pylith::topology::Field::Description description;
    description.label = subfieldName;
    description.alias = subfieldName;
    description.vectorFieldType = pylith::topology::Field::SCALAR;
    description.numComponents = 1;
    description.componentNames.resize(1);
    description.componentNames[0] = subfieldName;
    description.scale = _normalizer->getLengthScale();
    description.validator = NULL;

    _field->subfieldAdd(description, getSubfieldDiscretization(subfieldName));
    this->setSubfieldQuery(subfieldName);

    PYLITH_METHOD_END;
} // addWellboreRadius


// ----------------------------------------------------------------------------
// Add wellbore length subfield to auxiliary fields.
void
pylith::sources::AuxiliaryFactoryWellboreSource::addWellboreLength(void) { // wellboreLength
    PYLITH_METHOD_BEGIN;
    PYLITH_JOURNAL_DEBUG("addWellboreLength(void)");

    const char* subfieldName = "wellbore_length";

    const PylithReal lengthScale = _normalizer->getLengthScale();

    pylith::topology::Field::Description description;
    description.label = subfieldName;
    description.alias = subfieldName;
    description.vectorFieldType = pylith::topology::Field::SCALAR;
    description.numComponents = 1;
    description.componentNames.resize(1);
    description.componentNames[0] = subfieldName;
    description.scale = _normalizer->getLengthScale();
    description.validator = NULL;

    _field->subfieldAdd(description, getSubfieldDiscretization(subfieldName));
    this->setSubfieldQuery(subfieldName);

    PYLITH_METHOD_END;
} // addWellboreLength


// --------------------------------------------------------------------
// Add wellbore pressure subfield to auxiliary fields.
void
pylith::sources::AuxiliaryFactoryWellboreSource::addWellborePressure(void) { // wellborePressure
    PYLITH_METHOD_BEGIN;
    PYLITH_JOURNAL_DEBUG("addWellborePressure(void)");

    const char* subfieldName = "wellbore_pressure";
    const PylithReal pressureScale = _normalizer->getPressureScale();

    pylith::topology::Field::Description description;
    description.label = subfieldName;
    description.alias = subfieldName;
    description.vectorFieldType = pylith::topology::Field::SCALAR;
    description.numComponents = 1;
    description.componentNames.resize(1);
    description.componentNames[0] = subfieldName;
    description.scale = pressureScale;
    description.validator = NULL;

    _field->subfieldAdd(description, getSubfieldDiscretization(subfieldName));
    this->setSubfieldQuery(subfieldName);

    PYLITH_METHOD_END;
} // addWellborePressure


// ----------------------------------------------------------------------
// Add wellbore character subfield to auxiliary fields.
void
pylith::sources::AuxiliaryFactoryWellboreSource::addWellboreCharacter(void) { // wellboreCharacter
    PYLITH_METHOD_BEGIN;
    PYLITH_JOURNAL_DEBUG("addWellboreCharacter(void)");

    const char* subfieldName = "wellbore_character";
    const PylithReal noScale = 1;

    pylith::topology::Field::Description description;
    description.label = subfieldName;
    description.alias = subfieldName;
    description.vectorFieldType = pylith::topology::Field::SCALAR;
    description.numComponents = 1;
    description.componentNames.resize(1);
    description.componentNames[0] = subfieldName;
    description.scale = noScale;
    description.validator = NULL;

    _field->subfieldAdd(description, getSubfieldDiscretization(subfieldName));
    this->setSubfieldQuery(subfieldName);

    PYLITH_METHOD_END;
} // addWellboreCharacter


// ----------------------------------------------------------------------------
// Add isotropic permeability subfield to auxiliary fields.
void
pylith::sources::AuxiliaryFactoryWellboreSource::addElementDimensions(void) { // elementLength
    PYLITH_METHOD_BEGIN;
    PYLITH_JOURNAL_DEBUG("addElementDimensions(void)");

    const char* subfieldName = "element_dimensions";
    const char* componentNames[3] = {
        "element_x",
        "element_y",
        "element_z"
    };
    const int tensorSize = (3 == _spaceDim) ? 3 : (2 == _spaceDim) ? 2 : 1;
    const PylithReal lengthScale = _normalizer->getLengthScale();

    pylith::topology::Field::Description description;
    description.label = subfieldName;
    description.alias = subfieldName;
    description.vectorFieldType = pylith::topology::Field::OTHER;
    description.numComponents = tensorSize;
    description.componentNames.resize(tensorSize);
    for (int i = 0; i < _spaceDim; ++i) {
        description.componentNames[i] = componentNames[i];
    } // for
    description.scale = _normalizer->getLengthScale();
    description.validator = NULL;

    _field->subfieldAdd(description, getSubfieldDiscretization(subfieldName));
    this->setSubfieldQuery(subfieldName);

    PYLITH_METHOD_END;
} // addElementLength


// End of file
