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

#include "pylith/bc/DiagnosticFieldFactory.hh" // implementation of object methods

#include "pylith/topology/Field.hh" // USES Field

#include "spatialdata/units/Nondimensional.hh" // USES Nondimensional

#include "pylith/utils/error.hh" // USES PYLITH_METHOD*
#include "pylith/utils/journals.hh" // USES PYLITH_JOURNAL*

#include <cassert>

// ------------------------------------------------------------------------------------------------
// Default constructor.
pylith::bc::DiagnosticFieldFactory::DiagnosticFieldFactory(void) {
    GenericComponent::setName("diagnosticfieldfactory");
} // constructor


// ------------------------------------------------------------------------------------------------
// Destructor.
pylith::bc::DiagnosticFieldFactory::~DiagnosticFieldFactory(void) {}


// ------------------------------------------------------------------------------------------------
// Add fault normal direction subfield to derived fields.
void
pylith::bc::DiagnosticFieldFactory::addNormalDir(void) {
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
// Add (horizontal) tangential direction subfield to derived fields.
void
pylith::bc::DiagnosticFieldFactory::addTangentialDirHoriz(void) {
    PYLITH_METHOD_BEGIN;
    PYLITH_JOURNAL_DEBUG("addTangentialDirHoriz(void)");

    const char* fieldName = _spaceDim == 2 ? "tangential_dir" : "horizontal_tangential_dir";
    const char* componentNames2D[2] = {
        "tangential_dir_x",
        "tangential_dir_y",
    };
    const char* componentNames3D[3] = {
        "horizontal_tangential_dir_x",
        "horizontal_tangential_dir_y",
        "horizontal_tangential_dir_z",
    };

    pylith::topology::Field::Description description;
    description.label = fieldName;
    description.alias = fieldName;
    description.vectorFieldType = pylith::topology::Field::VECTOR;
    description.numComponents = _spaceDim;
    description.componentNames.resize(_spaceDim);
    for (int i = 0; i < _spaceDim; ++i) {
        description.componentNames[i] = _spaceDim == 2 ? componentNames2D[i] : componentNames3D[i];
    } // for
    description.scale = 1.0;
    description.validator = NULL;

    _field->subfieldAdd(description, getSubfieldDiscretization(fieldName));
    // No query; computed during initialization.

    PYLITH_METHOD_END;
} // addTangentialDirHoriz


// ------------------------------------------------------------------------------------------------
// Add (vertical) tangential direction subfield to derived fields.
void
pylith::bc::DiagnosticFieldFactory::addTangentialDirVert(void) {
    PYLITH_METHOD_BEGIN;
    PYLITH_JOURNAL_DEBUG("addTangentialDirVert(void)");

    const char* fieldName = "vertical_tangential_dir";
    const char* componentNames[3] = {
        "vertical_tangential_dir_x",
        "vertical_tangential_dir_y",
        "vertical_tangential_dir_z"
    };

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
} // addTangentialDirVert


// ------------------------------------------------------------------------------------------------
// Add subfields using discretizations provided.
void
pylith::bc::DiagnosticFieldFactory::addSubfields(void) {
    PYLITH_METHOD_BEGIN;
    PYLITH_JOURNAL_DEBUG("addSubfields(void)");

    if (_subfieldDiscretizations.find("normal_dir") != _subfieldDiscretizations.end()) {
        addNormalDir();
    } // if
    if ((_subfieldDiscretizations.find("tangential_dir") != _subfieldDiscretizations.end()) ||
        (_subfieldDiscretizations.find("horizontal_tangential_dir") != _subfieldDiscretizations.end())) {
        addTangentialDirHoriz();
    } // if
    if ((_spaceDim > 2) && (_subfieldDiscretizations.find("vertical_tangential_dir") != _subfieldDiscretizations.end())) {
        addTangentialDirVert();
    } // if

    PYLITH_METHOD_END;

}


// End of file
