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

#include "TimeDependentAuxiliaryFactory.hh" // implementation of object methods

#include "pylith/topology/Field.hh" // HOLDSA AuxiliaryField
#include "pylith/topology/FieldQuery.hh" // USES FieldQuery

#include "spatialdata/units/Nondimensional.hh" // USES Nondimensional

#include "pylith/utils/journals.hh" // USES PYLITH_JOURNAL*

#include <cassert>

// ----------------------------------------------------------------------
const char* pylith::bc::TimeDependentAuxiliaryFactory::_genericComponent = "timedependentauxiliaryfactory";

const char* pylith::bc::TimeDependentAuxiliaryFactory::_componentsXYZ[3] = { "_x", "_y", "_z" };
const char* pylith::bc::TimeDependentAuxiliaryFactory::_componentsTN[2] = { "_tangential", "_normal" };
const char* pylith::bc::TimeDependentAuxiliaryFactory::_componentsTTN[3] = { "_tangential_1", "_tangential_2", "_normal" };

// ----------------------------------------------------------------------
// Default constructor.
pylith::bc::TimeDependentAuxiliaryFactory::TimeDependentAuxiliaryFactory(const ReferenceEnum reference) :
    _auxComponents(reference)
{ // constructor
    GenericComponent::name(_genericComponent);
} // constructor

// ----------------------------------------------------------------------
// Destructor.
pylith::bc::TimeDependentAuxiliaryFactory::~TimeDependentAuxiliaryFactory(void)
{ // destructor
} // destructor

// ----------------------------------------------------------------------
// Set names of vector components in auxiliary subfield.
void
pylith::bc::TimeDependentAuxiliaryFactory::_setVectorFieldComponentNames(pylith::topology::FieldBase::Description* description) {
    PYLITH_METHOD_BEGIN;
    PYLITH_JOURNAL_DEBUG("_setVectorFieldComponentNames(description)");

    assert(description);

    const char** componentNames = NULL;
    if (XYZ == _auxComponents) {
        componentNames = _componentsXYZ;
    } else if (TANGENTIAL_NORMAL == _auxComponents) {
        if (2 == _spaceDim) {
            componentNames = _componentsTN;
        } else if (3 == _spaceDim) {
            componentNames = _componentsTTN;
        } // if/else
    } // if/else
    if (!componentNames) {
        PYLITH_JOURNAL_ERROR("Unknown case for auxiliary component reference ("<<_auxComponents<<") and spatial dimension ("<<_spaceDim<<").");
        throw std::logic_error("Unknown case for auxiliary component reference and spatial dimension.");
    } // if

    assert(size_t(_spaceDim) == description->numComponents);
    for (int i = 0; i < _spaceDim; ++i) {
        description->componentNames[i] = description->label + std::string(componentNames[i]);
    } // for

    PYLITH_METHOD_END;
} // _setVectorFieldComponentNames


// ----------------------------------------------------------------------
// Add initial amplitude field to auxiliary fields.
void
pylith::bc::TimeDependentAuxiliaryFactory::initialAmplitude(void)
{ // initialAmplitudert
    PYLITH_METHOD_BEGIN;
    PYLITH_JOURNAL_DEBUG("initialAmplitude(void)");

    const char* fieldName = "initial_amplitude";

    assert(_defaultDescription);
    pylith::topology::FieldBase::Description subfieldDescription(*_defaultDescription);

    subfieldDescription.label = fieldName;
    subfieldDescription.alias = fieldName;
    subfieldDescription.validator = NULL;
    switch (subfieldDescription.vectorFieldType) {
    case pylith::topology::FieldBase::SCALAR: {
        const size_t numComponents = 1;
        assert(numComponents == subfieldDescription.numComponents);
        assert(numComponents == subfieldDescription.componentNames.size());
        subfieldDescription.componentNames[0] = fieldName;
        break;
    } // SCALAR
    case pylith::topology::FieldBase::VECTOR: {
        _setVectorFieldComponentNames(&subfieldDescription);
        break;
    } // VECTOR
    default:
        PYLITH_JOURNAL_ERROR("Unknown vector field case.");
        throw std::logic_error("Unknown vector field case in TimeDependentAuxiliaryFactory::initialAmplitude().");
    } // switch

    _field->subfieldAdd(subfieldDescription, subfieldDiscretization(fieldName));
    _subfieldQueryFn(fieldName, pylith::topology::FieldQuery::dbQueryGeneric);

    PYLITH_METHOD_END;
} // initialAmplitude


// ----------------------------------------------------------------------
// Add rate amplitude field to auxiliary fields.
void
pylith::bc::TimeDependentAuxiliaryFactory::rateAmplitude(void)
{ // rateAmplitude
    PYLITH_METHOD_BEGIN;
    PYLITH_JOURNAL_DEBUG("rateAmplitude(void)");

    const char* fieldName = "rate_amplitude";

    assert(_defaultDescription);
    pylith::topology::FieldBase::Description subfieldDescription(*_defaultDescription);
    subfieldDescription.label = fieldName;
    subfieldDescription.alias = fieldName;
    subfieldDescription.validator = NULL;
    switch (subfieldDescription.vectorFieldType) {
    case pylith::topology::FieldBase::SCALAR: {
        const size_t numComponents = 1;
        assert(numComponents == subfieldDescription.numComponents);
        assert(numComponents == subfieldDescription.componentNames.size());
        subfieldDescription.componentNames[0] = fieldName;
        break;
    } // SCALAR
    case pylith::topology::FieldBase::VECTOR: {
        _setVectorFieldComponentNames(&subfieldDescription);
        break;
    } // VECTOR
    default:
        PYLITH_JOURNAL_ERROR("Unknown vector field case.");
        throw std::logic_error("Unknown vector field case in TimeDependentAuxiliaryFactory::rateAmplitude().");
    } // switch

    _field->subfieldAdd(subfieldDescription, subfieldDiscretization(fieldName));
    _subfieldQueryFn(fieldName, pylith::topology::FieldQuery::dbQueryGeneric);

    PYLITH_METHOD_END;
} // rateAmplitude


// ----------------------------------------------------------------------
// Add rate start time field to auxiliary fields.
void
pylith::bc::TimeDependentAuxiliaryFactory::rateStartTime(void)
{ // rateStartTime
    PYLITH_METHOD_BEGIN;
    PYLITH_JOURNAL_DEBUG("rateStartTime(void)");

    const char* fieldName = "rate_start_time";

    assert(_defaultDescription);
    assert(_normalizer);
    pylith::topology::FieldBase::Description subfieldDescription(*_defaultDescription);
    subfieldDescription.label = fieldName;
    subfieldDescription.alias = fieldName;
    subfieldDescription.vectorFieldType = pylith::topology::FieldBase::SCALAR;
    subfieldDescription.numComponents = 1;
    subfieldDescription.componentNames.resize(1);
    subfieldDescription.componentNames[0] = fieldName;
    subfieldDescription.scale = _normalizer->timeScale();
    subfieldDescription.validator = NULL;

    _field->subfieldAdd(subfieldDescription, subfieldDiscretization(fieldName));
    _subfieldQueryFn(fieldName, pylith::topology::FieldQuery::dbQueryGeneric);

    PYLITH_METHOD_END;
} // rateStartTime


// ----------------------------------------------------------------------
// Add time history amplitude field to auxiliary fields.
void
pylith::bc::TimeDependentAuxiliaryFactory::timeHistoryAmplitude(void)
{ // timeHistoryAmplitude
    PYLITH_METHOD_BEGIN;
    PYLITH_JOURNAL_DEBUG("timeHistoryAmplitude(void)");

    const char* fieldName = "time_history_amplitude";

    assert(_defaultDescription);
    pylith::topology::FieldBase::Description subfieldDescription(*_defaultDescription);

    subfieldDescription.label = fieldName;
    subfieldDescription.alias = fieldName;
    switch (subfieldDescription.vectorFieldType) {
    case pylith::topology::FieldBase::SCALAR: {
        const size_t numComponents = 1;
        assert(numComponents == subfieldDescription.numComponents);
        assert(numComponents == subfieldDescription.componentNames.size());
        subfieldDescription.componentNames[0] = fieldName;
        break;
    } // SCALAR
    case pylith::topology::FieldBase::VECTOR: {
        _setVectorFieldComponentNames(&subfieldDescription);
        break;
    } // VECTOR
    default:
        PYLITH_JOURNAL_ERROR("Unknown vector field case.");
        throw std::logic_error("Unknown vector field case in TimeDependentAuxiliaryFactory::timeHistoryAmplitude().");
    } // switch

    _field->subfieldAdd(subfieldDescription, subfieldDiscretization(fieldName));
    _subfieldQueryFn(fieldName, pylith::topology::FieldQuery::dbQueryGeneric);

    PYLITH_METHOD_END;
} // timeHistoryAmplitude


// ----------------------------------------------------------------------
// Add time history start time field to auxiliary fields.
void
pylith::bc::TimeDependentAuxiliaryFactory::timeHistoryStartTime(void)
{ // timeHistoryStartTime
    PYLITH_METHOD_BEGIN;
    PYLITH_JOURNAL_DEBUG("timeHistoryStartTime(void)");

    const char* fieldName = "time_history_start_time";

    assert(_defaultDescription);
    assert(_normalizer);
    pylith::topology::FieldBase::Description subfieldDescription(*_defaultDescription);
    subfieldDescription.label = fieldName;
    subfieldDescription.alias = fieldName;
    subfieldDescription.vectorFieldType = pylith::topology::FieldBase::SCALAR;
    subfieldDescription.numComponents = 1;
    subfieldDescription.componentNames.resize(1);
    subfieldDescription.componentNames[0] = fieldName;
    subfieldDescription.scale = _normalizer->timeScale();
    subfieldDescription.validator = NULL;

    _field->subfieldAdd(subfieldDescription, subfieldDiscretization(fieldName));
    _subfieldQueryFn(fieldName, pylith::topology::FieldQuery::dbQueryGeneric);

    PYLITH_METHOD_END;
} // timeHistoryStartTime


// ----------------------------------------------------------------------
// Add time history value field to auxiliary fields.
void
pylith::bc::TimeDependentAuxiliaryFactory::timeHistoryValue(void)
{ // timeHistoryValue
    PYLITH_METHOD_BEGIN;
    PYLITH_JOURNAL_DEBUG("timeHistoryValue(void)");

    const char* fieldName = "time_history_value";

    assert(_defaultDescription);
    assert(_normalizer);
    pylith::topology::FieldBase::Description subfieldDescription(*_defaultDescription);
    subfieldDescription.label = fieldName;
    subfieldDescription.alias = fieldName;
    subfieldDescription.vectorFieldType = pylith::topology::FieldBase::SCALAR;
    subfieldDescription.numComponents = 1;
    subfieldDescription.componentNames.resize(1);
    subfieldDescription.componentNames[0] = fieldName;
    subfieldDescription.scale = _normalizer->timeScale();
    subfieldDescription.validator = NULL;

    _field->subfieldAdd(subfieldDescription, subfieldDiscretization(fieldName));
    _subfieldQueryFn(fieldName, NULL); // populated by integrator or constraint at begining of time step.

    PYLITH_METHOD_END;
} // timeHistoryValue

// End of file
