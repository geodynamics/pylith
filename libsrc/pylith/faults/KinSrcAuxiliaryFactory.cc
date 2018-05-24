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

#include "KinSrcAuxiliaryFactory.hh" // implementation of object methods

#include "pylith/topology/Field.hh" // USES Field
#include "pylith/topology/FieldQuery.hh" // HOLDSA FieldQuery

#include "spatialdata/units/Nondimensional.hh" // USES Nondimensional

#include "pylith/utils/journals.hh" // USES PYLITH_JOURNAL*

#include <cassert>

// ----------------------------------------------------------------------
const char* pylith::faults::KinSrcAuxiliaryFactory::_genericComponent = "kinsrcauxiliaryfactory";

// ----------------------------------------------------------------------
// Default constructor.
pylith::faults::KinSrcAuxiliaryFactory::KinSrcAuxiliaryFactory(void)
{ // constructor
    GenericComponent::name(_genericComponent);
} // constructor

// ----------------------------------------------------------------------
// Destructor.
pylith::faults::KinSrcAuxiliaryFactory::~KinSrcAuxiliaryFactory(void) {}

// ----------------------------------------------------------------------
// Add slipTime subfield to auxiliary fields.
void
pylith::faults::KinSrcAuxiliaryFactory::slipTime(void) {
    PYLITH_METHOD_BEGIN;
    PYLITH_JOURNAL_DEBUG("slipTime(void)");

    const char* fieldName = "slip_time";
    const PylithReal timeScale = _normalizer->timeScale();

    pylith::topology::Field::Description description;
    description.label = fieldName;
    description.vectorFieldType = pylith::topology::Field::SCALAR;
    description.numComponents = 1;
    description.componentNames.resize(1);
    description.componentNames[0] = fieldName;
    description.scale = timeScale;
    description.validator = NULL;

    _field->subfieldAdd(description, subfieldDiscretization(fieldName));
    _subfieldQueryFn(fieldName, pylith::topology::FieldQuery::dbQueryGeneric);

    PYLITH_METHOD_END;
} // slipTime


// ----------------------------------------------------------------------
// Add riseTime subfield to auxiliary fields.
void
pylith::faults::KinSrcAuxiliaryFactory::riseTime(void) {
    PYLITH_METHOD_BEGIN;
    PYLITH_JOURNAL_DEBUG("riseTime(void)");

    const char* fieldName = "rise_time";
    const PylithReal timeScale = _normalizer->timeScale();

    pylith::topology::Field::Description description;
    description.label = fieldName;
    description.vectorFieldType = pylith::topology::Field::SCALAR;
    description.numComponents = 1;
    description.componentNames.resize(1);
    description.componentNames[0] = fieldName;
    description.scale = timeScale;
    description.validator = NULL;

    _field->subfieldAdd(description, subfieldDiscretization(fieldName));
    _subfieldQueryFn(fieldName, pylith::topology::FieldQuery::dbQueryGeneric);

    PYLITH_METHOD_END;
} // riseTime


// ----------------------------------------------------------------------
// Add final slip subfield to auxiliary fields.
void
pylith::faults::KinSrcAuxiliaryFactory::finalSlip(void) {
    PYLITH_METHOD_BEGIN;
    PYLITH_JOURNAL_DEBUG("finalSlip(void)");

    const char* fieldName = "final_slip";
    const char* componentNames2D[2] = { "final_slip_left_lateral", "final_slip_opening" };
    const char* componentNames3D[3] = { "final_slip_left_lateral", "final_slip_reverse", "final_slip_opening" };


    const PylithReal lengthScale = _normalizer->lengthScale();

    pylith::topology::Field::Description description;
    description.label = fieldName;
    description.vectorFieldType = pylith::topology::Field::VECTOR;
    description.numComponents = _spaceDim;
    description.componentNames.resize(_spaceDim);
    switch (_spaceDim) {
    case 2:
        for (int i = 0; i < 2; ++i) {
            description.componentNames[i] = componentNames2D[i];
        } // for
        break;
    case 3:
        for (int i = 0; i < 3; ++i) {
            description.componentNames[i] = componentNames3D[i];
        } // for
        break;
    default:
        PYLITH_JOURNAL_ERROR("Unknown spatial dimension ("<<_spaceDim<<") in setting up final_slip auxiliary subfield.");
        throw std::logic_error("Unknown spatial dimension in setting up final_slip auxiliary subfield.");
    } // switch
    description.scale = lengthScale;
    description.validator = NULL;

    _field->subfieldAdd(description, subfieldDiscretization(fieldName));
    _subfieldQueryFn(fieldName, pylith::topology::FieldQuery::dbQueryGeneric);

    PYLITH_METHOD_END;
} // finalSlip


// ----------------------------------------------------------------------
// Add slip rate subfield to auxiliary fields.
void
pylith::faults::KinSrcAuxiliaryFactory::slipRate(void) {
    PYLITH_METHOD_BEGIN;
    PYLITH_JOURNAL_DEBUG("slipRate(void)");

    const char* fieldName = "slip_rate";
    const char* componentNames2D[2] = { "slip_rate_left_lateral", "slip_rate_opening" };
    const char* componentNames3D[3] = { "slip_rate_left_lateral", "slip_rate_reverse", "slip_rate_opening" };

    const PylithReal lengthScale = _normalizer->lengthScale();
    const PylithReal timeScale = _normalizer->timeScale();

    pylith::topology::Field::Description description;
    description.label = fieldName;
    description.vectorFieldType = pylith::topology::Field::VECTOR;
    description.numComponents = _spaceDim;
    description.componentNames.resize(_spaceDim);
    switch (_spaceDim) {
    case 2:
        for (int i = 0; i < 2; ++i) {
            description.componentNames[i] = componentNames2D[i];
        } // for
        break;
    case 3:
        for (int i = 0; i < 3; ++i) {
            description.componentNames[i] = componentNames3D[i];
        } // for
        break;
    default:
        PYLITH_JOURNAL_ERROR("Unknown spatial dimension ("<<_spaceDim<<") in setting up slip_rate auxiliary subfield.");
        throw std::logic_error("Unknown spatial dimension in setting up slip_rate auxiliary subfield.");
    } // switch
    description.scale = lengthScale / timeScale;
    description.validator = NULL;

    _field->subfieldAdd(description, subfieldDiscretization(fieldName));
    _subfieldQueryFn(fieldName, pylith::topology::FieldQuery::dbQueryGeneric);

    PYLITH_METHOD_END;
} // slipRate


// End of file
