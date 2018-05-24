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
// Copyright (c) 2010-2017 University of California, Davis
//
// See COPYING for license information.
//
// ----------------------------------------------------------------------
//

#include <portinfo>

#include "AbsorbingDampersAuxiliaryFactory.hh" // implementation of object methods

#include "pylith/topology/Field.hh" // HOLDSA AuxiliaryField
#include "pylith/topology/FieldQuery.hh" // USES FieldQuery

#include "spatialdata/units/Nondimensional.hh" // USES Nondimensional

#include "pylith/utils/journals.hh" // USES PYLITH_JOURNAL*

#include <cassert>

// ----------------------------------------------------------------------
const char* pylith::bc::AbsorbingDampersAuxiliaryFactory::_genericComponent = "absorbingdampersauxiliaryfactory";

// ----------------------------------------------------------------------
// Default constructor.
pylith::bc::AbsorbingDampersAuxiliaryFactory::AbsorbingDampersAuxiliaryFactory(void)
{ // constructor
    GenericComponent::name(_genericComponent);
} // constructor

// ----------------------------------------------------------------------
// Destructor.
pylith::bc::AbsorbingDampersAuxiliaryFactory::~AbsorbingDampersAuxiliaryFactory(void) {}

// ----------------------------------------------------------------------
// Add density field to auxiliary fields.
void
pylith::bc::AbsorbingDampersAuxiliaryFactory::density(void) {
    PYLITH_METHOD_BEGIN;
    PYLITH_JOURNAL_DEBUG("density(void)");

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

    _field->subfieldAdd(description, subfieldDiscretization(fieldName));
    _subfieldQueryFn(fieldName, pylith::topology::FieldQuery::dbQueryGeneric);

    PYLITH_METHOD_END;
} // density


// ----------------------------------------------------------------------
// Add shear wave speed field to auxiliary fields.
void
pylith::bc::AbsorbingDampersAuxiliaryFactory::vs(void) {
    PYLITH_METHOD_BEGIN;
    PYLITH_JOURNAL_DEBUG("vs(void)");

    const char* fieldName = "vs";
    const PylithReal velocityScale = _normalizer->lengthScale() / _normalizer->timeScale();

    pylith::topology::Field::Description description;
    description.label = fieldName;
    description.alias = fieldName;
    description.vectorFieldType = pylith::topology::Field::SCALAR;
    description.numComponents = 1;
    description.componentNames.resize(1);
    description.componentNames[0] = fieldName;
    description.scale = velocityScale;
    description.validator = NULL;

    _field->subfieldAdd(description, subfieldDiscretization(fieldName));
    _subfieldQueryFn(fieldName, pylith::topology::FieldQuery::dbQueryGeneric);

    PYLITH_METHOD_END;
} // vs


// ----------------------------------------------------------------------
// Add dilatational wave speed field to auxiliary fields.
void
pylith::bc::AbsorbingDampersAuxiliaryFactory::vp(void) {
    PYLITH_METHOD_BEGIN;
    PYLITH_JOURNAL_DEBUG("vp(void)");

    const char* fieldName = "vp";
    const PylithReal velocityScale = _normalizer->lengthScale() / _normalizer->timeScale();

    pylith::topology::Field::Description description;
    description.label = fieldName;
    description.alias = fieldName;
    description.vectorFieldType = pylith::topology::Field::SCALAR;
    description.numComponents = 1;
    description.componentNames.resize(1);
    description.componentNames[0] = fieldName;
    description.scale = velocityScale;
    description.validator = NULL;

    _field->subfieldAdd(description, subfieldDiscretization(fieldName));
    _subfieldQueryFn(fieldName, pylith::topology::FieldQuery::dbQueryGeneric);

    PYLITH_METHOD_END;
} // vp


// End of file
