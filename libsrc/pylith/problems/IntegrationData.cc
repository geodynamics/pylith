// -*- C++ -*-
//
// ======================================================================
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
// See LICENSE.md for license information.
//
// ======================================================================
//

#include <portinfo>

#include "IntegrationData.hh" // implementation of class methods

#include "pylith/topology/Field.hh" // USES Field
#include "pylith/utils/journals.hh" // USES PYLITH_JOURNAL_*
#include "pylith/utils/error.hh" // USES PYLITH_METHOD_*

#include <cassert> // USES assert()
#include <typeinfo> // USES typeid()
#include <stdexcept> // USES std::runtime_error

// ------------------------------------------------------------------------------------------------
const std::string pylith::problems::IntegrationData::time = "t";
const std::string pylith::problems::IntegrationData::time_step = "dt";
const std::string pylith::problems::IntegrationData::s_tshift = "s_tshift";
const std::string pylith::problems::IntegrationData::t_state = "t_state";
const std::string pylith::problems::IntegrationData::dt_residual = "dt_residual";
const std::string pylith::problems::IntegrationData::dt_jacobian = "dt_jacobian";
const std::string pylith::problems::IntegrationData::dt_lumped_jacobian_inverse = "dt_lumped_jacobian_inverse";

const std::string pylith::problems::IntegrationData::solution = "solution";
const std::string pylith::problems::IntegrationData::solution_dot = "solution_dot";
const std::string pylith::problems::IntegrationData::residual = "residual";
const std::string pylith::problems::IntegrationData::lumped_jacobian_inverse = "lumped_jacobian_inverse";

// ------------------------------------------------------------------------------------------------
// Constructor
pylith::problems::IntegrationData::IntegrationData(void) {
    GenericComponent::setName("integrationdata");
}


// ------------------------------------------------------------------------------------------------
// Destructor
pylith::problems::IntegrationData::~IntegrationData(void) {
    deallocate();
} // destructor


// ------------------------------------------------------------------------------------------------

// Deallocate data.
void
pylith::problems::IntegrationData::deallocate(void) {
    _scalars.clear();

    for (fields_map_t::iterator iter = _fields.begin(); iter != _fields.end(); ++iter) {
        if (iter->first != solution) { // Solution memory management handled by Python.
            delete iter->second;iter->second = NULL;
        } // if
    } // for
    _fields.clear();
} // deallocate


// ------------------------------------------------------------------------------------------------
// Set scalar quantity.
void
pylith::problems::IntegrationData::setScalar(const std::string& name,
                                             const PylithReal value) {
    PYLITH_METHOD_BEGIN;
    PYLITH_JOURNAL_DEBUG("setScalar(name="<<name<<", value="<<value<<")");

    _scalars[name] = value;

    PYLITH_METHOD_END;
}


// ------------------------------------------------------------------------------------------------
// Get scalar quantity.
PylithReal
pylith::problems::IntegrationData::getScalar(const std::string& name) const {
    PYLITH_METHOD_BEGIN;
    PYLITH_JOURNAL_DEBUG("getScalar(name="<<name<<")");

    scalars_map_t::const_iterator iter = _scalars.find(name);
    if (iter == _scalars.end()) {
        PYLITH_JOURNAL_LOGICERROR("No scalar value '" << name << "' in integration data.");
    } // if

    PYLITH_METHOD_RETURN(iter->second);
}


// ------------------------------------------------------------------------------------------------
// Remove scalar quantity.
void
pylith::problems::IntegrationData::removeScalar(const std::string& name) {
    PYLITH_METHOD_BEGIN;
    PYLITH_JOURNAL_DEBUG("removeScalar(name="<<name<<")");

    if (_scalars.count(name)) {
        _scalars.erase(name);
    } // if

    PYLITH_METHOD_END;
}


// ------------------------------------------------------------------------------------------------
// Check if we have field with given name.
bool
pylith::problems::IntegrationData::hasField(const std::string& name) const {
    return _fields.count(name) > 0;
}


// ------------------------------------------------------------------------------------------------
// Set field.
void
pylith::problems::IntegrationData::setField(const std::string& name,
                                            pylith::topology::Field* const field) {
    PYLITH_METHOD_BEGIN;
    PYLITH_JOURNAL_DEBUG("setField(name="<<name<<", field="<<typeid(field).name()<<")");

    fields_map_t::iterator iter = _fields.find(name);
    if (iter != _fields.end()) {
        delete iter->second;iter->second = field;
    } else {
        _fields[name] = field;
    } // if/else

    PYLITH_METHOD_END;
}


// ------------------------------------------------------------------------------------------------
// Get field.
pylith::topology::Field*
pylith::problems::IntegrationData::getField(const std::string& name) const {
    PYLITH_METHOD_BEGIN;
    PYLITH_JOURNAL_DEBUG("getField(name="<<name<<")");

    fields_map_t::const_iterator iter = _fields.find(name);
    if (iter == _fields.end()) {
        PYLITH_JOURNAL_LOGICERROR("No field '" << name << "' in integration data.");
    } // if

    PYLITH_METHOD_RETURN(iter->second);
}


// ------------------------------------------------------------------------------------------------
// Dump integration data to std::ostream.
std::string
pylith::problems::IntegrationData::str(void) const {
    std::ostringstream info;
    info << "Scalars:";
    for (scalars_map_t::const_iterator iter = _scalars.begin(); iter != _scalars.end(); ++iter) {
        info << " " << iter->first << "=" << iter->second;
    } // for
    info << "; ";

    info << "Fields:";
    for (fields_map_t::const_iterator iter = _fields.begin(); iter != _fields.end(); ++iter) {
        info << " " << iter->first;
    } // for
    info << "\n";

    return info.str();
}


// End of file
