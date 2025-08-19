// =================================================================================================
// This code is part of PyLith, developed through the Computational Infrastructure
// for Geodynamics (https://github.com/geodynamics/pylith).
//
// Copyright (c) 2010-2025, University of California, Davis and the PyLith Development Team.
// All rights reserved.
//
// See https://mit-license.org/ and LICENSE.md and for license information.
// =================================================================================================

#include <portinfo>

#include "pylith/faults/DiagnosticFieldFactory.hh" // implementation of object methods

#include "pylith/topology/Field.hh" // USES Field

#include "spatialdata/units/Scales.hh" // USES Scales

#include "pylith/utils/error.hh" // USES PYLITH_METHOD*
#include "pylith/utils/journals.hh" // USES PYLITH_JOURNAL*

#include <cassert>

// ------------------------------------------------------------------------------------------------
// Default constructor.
pylith::faults::DiagnosticFieldFactory::DiagnosticFieldFactory(void) {
    GenericComponent::setName("diagnosticfieldfactory");
} // constructor


// ------------------------------------------------------------------------------------------------
// Destructor.
pylith::faults::DiagnosticFieldFactory::~DiagnosticFieldFactory(void) {}


// ------------------------------------------------------------------------------------------------
// Add fault normal direction subfield to derived fields.
void
pylith::faults::DiagnosticFieldFactory::addNormalDir(void) {
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
    // No query; computed during initialization.

    PYLITH_METHOD_END;
} // addNormalDir


// ------------------------------------------------------------------------------------------------
// Add fault strike direction subfield to derived fields.
void
pylith::faults::DiagnosticFieldFactory::addStrikeDir(void) {
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
    // No query; computed during initialization.

    PYLITH_METHOD_END;
} // addStrikeDir


// ------------------------------------------------------------------------------------------------
// Add fault up-dip direction subfield to derived fields.
void
pylith::faults::DiagnosticFieldFactory::addUpDipDir(void) {
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
    // No query; computed during initialization.

    PYLITH_METHOD_END;
} // addUpDipDir


// ------------------------------------------------------------------------------------------------
// Add subfields using discretizations provided.
void
pylith::faults::DiagnosticFieldFactory::addSubfields(void) {
    PYLITH_METHOD_BEGIN;
    PYLITH_JOURNAL_DEBUG("addSubfields(void)");

    if (_subfieldDiscretizations.find("normal_dir") != _subfieldDiscretizations.end()) {
        addNormalDir();
    } // if
    if (_subfieldDiscretizations.find("strike_dir") != _subfieldDiscretizations.end()) {
        addStrikeDir();
    } // if
    if ((_spaceDim > 2) && (_subfieldDiscretizations.find("up_dip_dir") != _subfieldDiscretizations.end())) {
        addUpDipDir();
    } // if

    PYLITH_METHOD_END;

}


// End of file
