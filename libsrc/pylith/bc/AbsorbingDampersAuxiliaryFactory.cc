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
// See LICENSE.md.md.md for license information.
//
// ----------------------------------------------------------------------
//

#include <portinfo>

#include "AbsorbingDampersAuxiliaryFactory.hh" // implementation of object methods

#include "pylith/topology/Field.hh" // HOLDSA AuxiliaryField
#include "pylith/topology/FieldQuery.hh" // USES FieldQuery

#include "spatialdata/units/Nondimensional.hh" // USES Nondimensional

#include "pylith/utils/error.hh" // USES PYLITH_METHOD*
#include "pylith/utils/journals.hh" // USES PYLITH_JOURNAL*

#include <cassert>

// ---------------------------------------------------------------------------------------------------------------------
// Default constructor.
pylith::bc::AbsorbingDampersAuxiliaryFactory::AbsorbingDampersAuxiliaryFactory(void) {
    GenericComponent::setName("absorbingdampersauxiliaryfactory");
} // constructor


// ---------------------------------------------------------------------------------------------------------------------
// Destructor.
pylith::bc::AbsorbingDampersAuxiliaryFactory::~AbsorbingDampersAuxiliaryFactory(void) {}


// ---------------------------------------------------------------------------------------------------------------------
// Add density field to auxiliary fields.
void
pylith::bc::AbsorbingDampersAuxiliaryFactory::addDensity(void) {
    PYLITH_METHOD_BEGIN;
    PYLITH_JOURNAL_DEBUG("addDensity(void)");

    const char* subfieldName = "density";
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
} // addDensity


// ---------------------------------------------------------------------------------------------------------------------
// Add shear wave speed field to auxiliary fields.
void
pylith::bc::AbsorbingDampersAuxiliaryFactory::addVs(void) {
    PYLITH_METHOD_BEGIN;
    PYLITH_JOURNAL_DEBUG("addVs(void)");

    const char* subfieldName = "vs";
    const PylithReal velocityScale = _normalizer->getLengthScale() / _normalizer->getTimeScale();

    pylith::topology::Field::Description description;
    description.label = subfieldName;
    description.alias = subfieldName;
    description.vectorFieldType = pylith::topology::Field::SCALAR;
    description.numComponents = 1;
    description.componentNames.resize(1);
    description.componentNames[0] = subfieldName;
    description.scale = velocityScale;
    description.validator = NULL;

    _field->subfieldAdd(description, getSubfieldDiscretization(subfieldName));
    this->setSubfieldQuery(subfieldName);

    PYLITH_METHOD_END;
} // addVs


// ---------------------------------------------------------------------------------------------------------------------
// Add dilatational wave speed field to auxiliary fields.
void
pylith::bc::AbsorbingDampersAuxiliaryFactory::addVp(void) {
    PYLITH_METHOD_BEGIN;
    PYLITH_JOURNAL_DEBUG("addVp(void)");

    const char* subfieldName = "vp";
    const PylithReal velocityScale = _normalizer->getLengthScale() / _normalizer->getTimeScale();

    pylith::topology::Field::Description description;
    description.label = subfieldName;
    description.alias = subfieldName;
    description.vectorFieldType = pylith::topology::Field::SCALAR;
    description.numComponents = 1;
    description.componentNames.resize(1);
    description.componentNames[0] = subfieldName;
    description.scale = velocityScale;
    description.validator = NULL;

    _field->subfieldAdd(description, getSubfieldDiscretization(subfieldName));
    this->setSubfieldQuery(subfieldName);

    PYLITH_METHOD_END;
} // addVp


// End of file
