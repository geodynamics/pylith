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

#include "AuxiliaryFactoryPointForce.hh" // implementation of object methods

#include "pylith/topology/Field.hh" // USES Field
#include "pylith/topology/FieldQuery.hh" // HOLDSA FieldQuery

#include "spatialdata/units/Nondimensional.hh" // USES Nondimensional

#include "pylith/utils/error.hh" // USES PYLITH_METHOD*
#include "pylith/utils/journals.hh" // USES PYLITH_JOURNAL*

#include <cassert>

// ---------------------------------------------------------------------------------------------------------------------
// Default constructor.
pylith::sources::AuxiliaryFactoryPointForce::AuxiliaryFactoryPointForce(void) {
    GenericComponent::setName("auxiliaryfactorypointforce");
} // constructor


// ---------------------------------------------------------------------------------------------------------------------
// Destructor.
pylith::sources::AuxiliaryFactoryPointForce::~AuxiliaryFactoryPointForce(void) {}


// ----------------------------------------------------------------------
// Add fluid density subfield to auxiliary fields.
void
pylith::sources::AuxiliaryFactoryPointForce::addPointForce(void) {
    PYLITH_METHOD_BEGIN;
    PYLITH_JOURNAL_DEBUG("addPointForce(void)");

    const char *subfieldName = "point_force";
    const char* componentNames[3] = {
        "point_force_x",
        "point_force_y",
        "point_force_z",
    };
    const PylithReal pressureScale = _normalizer->getPressureScale();

    pylith::topology::Field::Description description;
    description.label = subfieldName;
    description.alias = subfieldName;
    description.vectorFieldType = pylith::topology::Field::VECTOR;
    description.numComponents = _spaceDim;
    description.componentNames.resize(_spaceDim);
    for (int i = 0; i < _spaceDim; ++i) {
        description.componentNames[i] = componentNames[i];
    } // for
    description.scale = pressureScale;
    description.validator = pylith::topology::FieldQuery::validatorPositive;

    _field->subfieldAdd(description, getSubfieldDiscretization(subfieldName));
    this->setSubfieldQuery(subfieldName);

    PYLITH_METHOD_END;
} // addPointForce


// ---------------------------------------------------------------------------------------------------------------------
// Add time delay of source time function to auxiliary fields.
void
pylith::sources::AuxiliaryFactoryPointForce::addTimeDelay(void) { // timeDelay
    PYLITH_METHOD_BEGIN;
    PYLITH_JOURNAL_DEBUG("addTimeDelay(void)");

    const char* subfieldName = "time_delay";

    pylith::topology::Field::Description description;
    const PylithReal timeScale = _normalizer->getTimeScale();
    description.label = subfieldName;
    description.alias = subfieldName;
    description.vectorFieldType = pylith::topology::Field::SCALAR;
    description.numComponents = 1;
    description.componentNames.resize(1);
    description.componentNames[0] = subfieldName;
    description.scale = timeScale;
    description.validator = pylith::topology::FieldQuery::validatorNonnegative;

    _field->subfieldAdd(description, getSubfieldDiscretization(subfieldName));
    this->setSubfieldQuery(subfieldName);

    PYLITH_METHOD_END;
} // addTimeDelay


// End of file
