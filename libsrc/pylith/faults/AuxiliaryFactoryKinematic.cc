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

#include "AuxiliaryFactoryKinematic.hh" // implementation of object methods

#include "pylith/topology/Field.hh" // USES Field
#include "pylith/topology/FieldQuery.hh" // HOLDSA FieldQuery

#include "spatialdata/units/Nondimensional.hh" // USES Nondimensional

#include "pylith/utils/journals.hh" // USES PYLITH_JOURNAL*

#include <cassert>

// ---------------------------------------------------------------------------------------------------------------------
// Default constructor.
pylith::faults::AuxiliaryFactoryKinematic::AuxiliaryFactoryKinematic(void) {
    GenericComponent::setName("auxiliaryfactorykinematic");
} // constructor


// ---------------------------------------------------------------------------------------------------------------------
// Destructor.
pylith::faults::AuxiliaryFactoryKinematic::~AuxiliaryFactoryKinematic(void) {}


// ---------------------------------------------------------------------------------------------------------------------
// Add fault strike direction subfield to auxiliary fields.
void
pylith::faults::AuxiliaryFactoryKinematic::addStrikeDir(void) {
    PYLITH_METHOD_BEGIN;
    PYLITH_JOURNAL_DEBUG("addStrikeDir(void)");

    const char* fieldName = "strike_dir";
    const char* componentNames[3] = { "strike_dir_x", "strike_dir_y", "strike_dir_z" };

    pylith::topology::Field::Description description;
    description.label = fieldName;
    description.alias = fieldName;
    description.vectorFieldType = pylith::topology::Field::VECTOR;
    description.numComponents = _spaceDim;
    description.componentNames.resize(_spaceDim);
    for (int i = 0; i < _spaceDim; ++i) {
        description.componentNames[i] = componentNames[i];
    } // for
    description.scale = 1.0;
    description.validator = NULL;

    _field->subfieldAdd(description, getSubfieldDiscretization(fieldName));
    _setSubfieldQueryFn(fieldName, NULL); // Computed during initialization.

    PYLITH_METHOD_END;
} // addStrikeDir


// ---------------------------------------------------------------------------------------------------------------------
// Add fault up-dip direction subfield to auxiliary fields.
void
pylith::faults::AuxiliaryFactoryKinematic::addUpDipDir(void) {
    PYLITH_METHOD_BEGIN;
    PYLITH_JOURNAL_DEBUG("addUpDipDir(void)");

    const char* fieldName = "up_dip_dir";
    const char* componentNames[3] = { "up_dip_dir_x", "up_dip_dir_y", "up_dip_dir_z" };

    pylith::topology::Field::Description description;
    description.label = fieldName;
    description.alias = fieldName;
    description.vectorFieldType = pylith::topology::Field::VECTOR;
    description.numComponents = _spaceDim;
    description.componentNames.resize(_spaceDim);
    for (int i = 0; i < _spaceDim; ++i) {
        description.componentNames[i] = componentNames[i];
    } // for
    description.scale = 1.0;
    description.validator = NULL;

    _field->subfieldAdd(description, getSubfieldDiscretization(fieldName));
    _setSubfieldQueryFn(fieldName, NULL); // Computed during initialization.

    PYLITH_METHOD_END;
} // addUpDipDir


// ---------------------------------------------------------------------------------------------------------------------
// Add fault normal direction subfield to auxiliary fields.
void
pylith::faults::AuxiliaryFactoryKinematic::addNormalDir(void) {
    PYLITH_METHOD_BEGIN;
    PYLITH_JOURNAL_DEBUG("addNormalDir(void)");

    const char* fieldName = "normal_dir";
    const char* componentNames[3] = { "normal_dir_x", "normal_dir_y", "normal_dir_z" };

    pylith::topology::Field::Description description;
    description.label = fieldName;
    description.alias = fieldName;
    description.vectorFieldType = pylith::topology::Field::VECTOR;
    description.numComponents = _spaceDim;
    description.componentNames.resize(_spaceDim);
    for (int i = 0; i < _spaceDim; ++i) {
        description.componentNames[i] = componentNames[i];
    } // for
    description.scale = 1.0;
    description.validator = NULL;

    _field->subfieldAdd(description, getSubfieldDiscretization(fieldName));
    _setSubfieldQueryFn(fieldName, NULL); // Computed during initialization.

    PYLITH_METHOD_END;
} // addNormalDir


// ---------------------------------------------------------------------------------------------------------------------
// Add fault slip subfield to auxiliary fields.
void
pylith::faults::AuxiliaryFactoryKinematic::addSlip(void) {
    PYLITH_METHOD_BEGIN;
    PYLITH_JOURNAL_DEBUG("addSlip(void)");

    const char* fieldName = "slip";
    const char* componentNames[3] = { "slip_opening", "slip_left_lateral", "slip_reverse" };

    const PylithReal lengthScale = _normalizer->getLengthScale();

    pylith::topology::Field::Description description;
    description.label = fieldName;
    description.alias = fieldName;
    description.vectorFieldType = pylith::topology::Field::VECTOR;
    description.numComponents = _spaceDim;
    description.componentNames.resize(_spaceDim);
    for (int i = 0; i < _spaceDim; ++i) {
        description.componentNames[i] = componentNames[i];
    } // for
    description.scale = lengthScale;
    description.validator = NULL;

    _field->subfieldAdd(description, getSubfieldDiscretization(fieldName));
    _setSubfieldQueryFn(fieldName, NULL); // Populated by kinematic source at beginning of time step.

    PYLITH_METHOD_END;
} // addSlip


// End of file
