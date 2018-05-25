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

#include "pylith/topology/Field.hh" // USES Field
#include "pylith/topology/FieldQuery.hh" // HOLDSA FieldQuery

#include "spatialdata/units/Nondimensional.hh" // USES Nondimensional

#include "pylith/utils/journals.hh" // USES PYLITH_JOURNAL*

#include <cassert>

// ----------------------------------------------------------------------
const char* pylith::faults::AuxiliaryFactory::_genericComponent = "faultsauxiliaryfactory";

// ----------------------------------------------------------------------
// Default constructor.
pylith::faults::AuxiliaryFactory::AuxiliaryFactory(void)
{ // constructor
    GenericComponent::name(_genericComponent);
} // constructor

// ----------------------------------------------------------------------
// Destructor.
pylith::faults::AuxiliaryFactory::~AuxiliaryFactory(void) {}

// ----------------------------------------------------------------------
// Add fault strike direction subfield to auxiliary fields.
void
pylith::faults::AuxiliaryFactory::strikeDir(void) {
    PYLITH_METHOD_BEGIN;
    PYLITH_JOURNAL_DEBUG("strikeDir(void)");

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

    _field->subfieldAdd(description, subfieldDiscretization(fieldName));
    _subfieldQueryFn(fieldName, NULL); // Computed during initialization.

    PYLITH_METHOD_END;
} // strikeDir

// ----------------------------------------------------------------------
// Add fault up-dip direction subfield to auxiliary fields.
void
pylith::faults::AuxiliaryFactory::upDipDir(void) {
    PYLITH_METHOD_BEGIN;
    PYLITH_JOURNAL_DEBUG("upDipDir(void)");

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

    _field->subfieldAdd(description, subfieldDiscretization(fieldName));
    _subfieldQueryFn(fieldName, NULL); // Computed during initialization.

    PYLITH_METHOD_END;
} // upDipDir

// ----------------------------------------------------------------------
// Add fault normal direction subfield to auxiliary fields.
void
pylith::faults::AuxiliaryFactory::normalDir(void) {
    PYLITH_METHOD_BEGIN;
    PYLITH_JOURNAL_DEBUG("normalDir(void)");

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

    _field->subfieldAdd(description, subfieldDiscretization(fieldName));
    _subfieldQueryFn(fieldName, NULL); // Computed during initialization.

    PYLITH_METHOD_END;
} // normalDir

// ----------------------------------------------------------------------
// Add fault slip subfield to auxiliary fields.
void
pylith::faults::AuxiliaryFactory::slip(void) {
    PYLITH_METHOD_BEGIN;
    PYLITH_JOURNAL_DEBUG("slip(void)");

    const char* fieldName = "slip";
    const char* componentNames[3] = { "slip_opening", "slip_left_lateral", "slip_reverse" };

    const PylithReal lengthScale = _normalizer->lengthScale();

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

    _field->subfieldAdd(description, subfieldDiscretization(fieldName));
    _subfieldQueryFn(fieldName, NULL); // Populated by kinematic source at beginning of time step.

    PYLITH_METHOD_END;
} // slip


// End of file
