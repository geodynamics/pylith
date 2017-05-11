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

#include "DirichletAuxiliaryFactory.hh" // implementation of object methods

#include "DirichletNew.hh" // USES DirichletNew

#include "pylith/topology/Field.hh" // USES Field
#include "pylith/topology/FieldQuery.hh" // HOLDSA FieldQuery

#include "spatialdata/units/Nondimensional.hh" // USES Nondimensional

#include "pylith/utils/journals.hh" // USES PYLITH_JOURNAL*

#include <cassert>

// ----------------------------------------------------------------------
const char* pylith::bc::DirichletAuxiliaryFactory::_genericComponent = "auxiliaryfactory";

// ----------------------------------------------------------------------
// Default constructor.
pylith::bc::DirichletAuxiliaryFactory::DirichletAuxiliaryFactory(const DirichletNew& bc,
								 const pylith::topology::Field& solution,
								 const PylithReal timeScale) :
    _bc(bc),
    _description(solution.subfieldInfo(bc._field.c_str()).description),
    _spaceDim(solution.spaceDim()),
    _timeScale(timeScale)
{ // constructor
    assert(1 <= _spaceDim && _spaceDim <= 3);
    assert(timeScale > 0.0);
    GenericComponent::name(_genericComponent);
} // constructor

// ----------------------------------------------------------------------
// Destructor.
pylith::bc::DirichletAuxiliaryFactory::~DirichletAuxiliaryFactory(void)
{ // destructor
} // destructor

// ----------------------------------------------------------------------
// Add initial amplitude field to auxiliary fields.
void
pylith::bc::DirichletAuxiliaryFactory::initialAmplitude(void) const
{ // initialAmplitude
    PYLITH_METHOD_BEGIN;
    PYLITH_JOURNAL_DEBUG("initialAmplitude(void)");

    const char* fieldName = "initial_amplitude";

    pylith::topology::FieldBase::Description auxDescription(_description);
    
    auxDescription.label = fieldName;
    auxDescription.validator = NULL;
    switch(auxDescription.vectorFieldType) {
    case pylith::topology::FieldBase::SCALAR: {
	const size_t numComponents = 1;
	assert(numComponents == auxDescription.numComponents);
	assert(numComponents == auxDescription.componentNames.size());
	auxDescription.componentNames[0] = "initial_amplitude";
	break;
    } // SCALAR	
    case pylith::topology::FieldBase::VECTOR: {
	const char* componentNames[3] = { "initial_amplitude_x", "initial_amplitude_y", "initial_amplitude_z" };
	assert(size_t(_spaceDim) == auxDescription.numComponents);
	assert(size_t(_spaceDim) == auxDescription.componentNames.size());
	for (int i=0; i < _spaceDim; ++i) {
	    auxDescription.componentNames[i] = componentNames[i];
	} // for
	break;
    } // VECTOR
	default :
	    PYLITH_JOURNAL_ERROR("Unknown vector field case.");
	    throw std::logic_error("Unknown vector field case in DirichletAuxiliaryFactory::initialAmplitude().");
    } // switch
    
    assert(_bc._auxFields);
    assert(_bc._auxFieldsQuery);
    _bc._auxFields->subfieldAdd(auxDescription, _bc.auxFieldDiscretization(fieldName));
    _bc._auxFieldsQuery->queryFn(fieldName, pylith::topology::FieldQuery::dbQueryGeneric);

    PYLITH_METHOD_END;
} // initialAmplitude


// ----------------------------------------------------------------------
// Add rate amplitude field to auxiliary fields.
void
pylith::bc::DirichletAuxiliaryFactory::rateAmplitude(void) const
{ // rateAmplitude
    PYLITH_METHOD_BEGIN;
    PYLITH_JOURNAL_DEBUG("rateAmplitude(void)");

    const char* fieldName = "rate_amplitude";

    pylith::topology::FieldBase::Description auxDescription(_description);
    auxDescription.label = fieldName;
    auxDescription.validator = NULL;
    switch(auxDescription.vectorFieldType) {
    case pylith::topology::FieldBase::SCALAR: {
	const size_t numComponents = 1;
	assert(numComponents == auxDescription.numComponents);
	assert(numComponents == auxDescription.componentNames.size());
	auxDescription.componentNames[0] = "rate_amplitude";
	break;
    } // SCALAR	
    case pylith::topology::FieldBase::VECTOR: {
	const char* componentNames[3] = { "rate_amplitude_x", "rate_amplitude_y", "rate_amplitude_z" };
	assert(size_t(_spaceDim) == auxDescription.numComponents);
	assert(size_t(_spaceDim) == auxDescription.componentNames.size());
	for (int i=0; i < _spaceDim; ++i) {
	    auxDescription.componentNames[i] = componentNames[i];
	} // for
	break;
    } // VECTOR
	default :
	    PYLITH_JOURNAL_ERROR("Unknown vector field case.");
	    throw std::logic_error("Unknown vector field case in DirichletAuxiliaryFactory::rateAmplitude().");
    } // switch
    
    assert(_bc._auxFields);
    assert(_bc._auxFieldsQuery);
    _bc._auxFields->subfieldAdd(auxDescription, _bc.auxFieldDiscretization(fieldName));
    _bc._auxFieldsQuery->queryFn(fieldName, pylith::topology::FieldQuery::dbQueryGeneric);

    PYLITH_METHOD_END;
} // rateAmplitude


// ----------------------------------------------------------------------
// Add rate start time field to auxiliary fields.
void
pylith::bc::DirichletAuxiliaryFactory::rateStartTime(void) const
{ // rateStartTime
    PYLITH_METHOD_BEGIN;
    PYLITH_JOURNAL_DEBUG("rateStartTime(void)");

    const char* fieldName = "rate_start_time";

    pylith::topology::FieldBase::Description auxDescription(_description);
    auxDescription.label = fieldName;
    auxDescription.vectorFieldType = pylith::topology::FieldBase::SCALAR;
    auxDescription.numComponents = 1;
    auxDescription.componentNames.resize(1);
    auxDescription.componentNames[0] = "rate_start_time";
    auxDescription.scale = _timeScale;
    auxDescription.validator = NULL;
    
    assert(_bc._auxFields);
    assert(_bc._auxFieldsQuery);
    _bc._auxFields->subfieldAdd(auxDescription, _bc.auxFieldDiscretization(fieldName));
    _bc._auxFieldsQuery->queryFn(fieldName, pylith::topology::FieldQuery::dbQueryGeneric);

    PYLITH_METHOD_END;
} // rateStartTime


// ----------------------------------------------------------------------
// Add time history amplitude field to auxiliary fields.
void
pylith::bc::DirichletAuxiliaryFactory::timeHistoryAmplitude(void) const
{ // timeHistoryAmplitude
    PYLITH_METHOD_BEGIN;
    PYLITH_JOURNAL_DEBUG("timeHistoryAmplitude(void)");

    const char* fieldName = "time_history_amplitude";

    pylith::topology::FieldBase::Description auxDescription(_description);
    
    auxDescription.label = fieldName;
    switch(auxDescription.vectorFieldType) {
    case pylith::topology::FieldBase::SCALAR: {
	const size_t numComponents = 1;
	assert(numComponents == auxDescription.numComponents);
	assert(numComponents == auxDescription.componentNames.size());
	auxDescription.componentNames[0] = "time_history_amplitude";
	break;
    } // SCALAR	
    case pylith::topology::FieldBase::VECTOR: {
	const char* componentNames[3] = { "time_history_amplitude_x", "time_history_amplitude_y", "time_history_amplitude_z" };
	assert(size_t(_spaceDim) == auxDescription.numComponents);
	assert(size_t(_spaceDim) == auxDescription.componentNames.size());
	for (int i=0; i < _spaceDim; ++i) {
	    auxDescription.componentNames[i] = componentNames[i];
	} // for
	break;
    } // VECTOR
	default :
	    PYLITH_JOURNAL_ERROR("Unknown vector field case.");
	    throw std::logic_error("Unknown vector field case in DirichletAuxiliaryFactory::timeHistoryAmplitude().");
    } // switch
    
    assert(_bc._auxFields);
    assert(_bc._auxFieldsQuery);
    _bc._auxFields->subfieldAdd(auxDescription, _bc.auxFieldDiscretization(fieldName));
    _bc._auxFieldsQuery->queryFn(fieldName, pylith::topology::FieldQuery::dbQueryGeneric);

    PYLITH_METHOD_END;
} // timeHistoryAmplitude


// ----------------------------------------------------------------------
// Add time history start time field to auxiliary fields.
void
pylith::bc::DirichletAuxiliaryFactory::timeHistoryStartTime(void) const
{ // timeHistoryStartTime
    PYLITH_METHOD_BEGIN;
    PYLITH_JOURNAL_DEBUG("timeHistoryStartTime(void)");

    const char* fieldName = "time_history_start_time";

    pylith::topology::FieldBase::Description auxDescription(_description);
    auxDescription.label = fieldName;
    auxDescription.vectorFieldType = pylith::topology::FieldBase::SCALAR;
    auxDescription.numComponents = 1;
    auxDescription.componentNames.resize(1);
    auxDescription.componentNames[0] = "time_history_start_time";
    auxDescription.scale = _timeScale;
    auxDescription.validator = NULL;
    
    assert(_bc._auxFields);
    assert(_bc._auxFieldsQuery);
    _bc._auxFields->subfieldAdd(auxDescription, _bc.auxFieldDiscretization(fieldName));
    _bc._auxFieldsQuery->queryFn(fieldName, pylith::topology::FieldQuery::dbQueryGeneric);

    PYLITH_METHOD_END;
} // timeHistoryStartTime


// ----------------------------------------------------------------------
// Add time history value field to auxiliary fields.
void
pylith::bc::DirichletAuxiliaryFactory::timeHistoryValue(void) const
{ // timeHistoryValue
    PYLITH_METHOD_BEGIN;
    PYLITH_JOURNAL_DEBUG("timeHistoryValue(void)");

    const char* fieldName = "time_history_value";

    pylith::topology::FieldBase::Description auxDescription(_description);
    auxDescription.label = fieldName;
    auxDescription.vectorFieldType = pylith::topology::FieldBase::SCALAR;
    auxDescription.numComponents = 1;
    auxDescription.componentNames.resize(1);
    auxDescription.componentNames[0] = "time_history_value";
    auxDescription.scale = _timeScale;
    auxDescription.validator = NULL;
    
    assert(_bc._auxFields);
    assert(_bc._auxFieldsQuery);
    _bc._auxFields->subfieldAdd(auxDescription, _bc.auxFieldDiscretization(fieldName));
    _bc._auxFieldsQuery->queryFn(fieldName, NULL); // populated by boundary condition

    PYLITH_METHOD_END;
} // timeHistoryValue


// End of file
