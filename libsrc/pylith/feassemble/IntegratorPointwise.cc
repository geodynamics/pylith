// -*- C++ -*-
//
// ======================================================================
//
// Brad T. Aagaard, U.S. Geological Survey
// Charles A. Williams, GNS Science
// Matthew G. Knepley, Rice University
//
// This code was developed as part of the Computational Infrastructure
// for Geodynamics (http://geodynamics.org).
//
// Copyright (c) 2010-2015 University of California, Davis
//
// See COPYING for license information.
//
// ======================================================================
//

#include <portinfo>

#include "IntegratorPointwise.hh" // implementation of class methods

#include "pylith/topology/Mesh.hh" // USES Mesh
#include "pylith/topology/Field.hh" // USES Field
#include "pylith/topology/FieldQuery.hh" // USES FieldQuery

#include "spatialdata/units/Nondimensional.hh" // USES Nondimensional

#include "pylith/utils/EventLogger.hh" // USES EventLogger

#include "journal/debug.h" // USES journal::debug_t

#include <cassert> // USES assert()
#include <stdexcept> // USES std::runtime_error

// ----------------------------------------------------------------------
// Constructor
pylith::feassemble::IntegratorPointwise::IntegratorPointwise(void) :
    _normalizer(0),
    _logger(0),
    _auxFields(0),
    _auxFieldsDB(0),
    _auxFieldsQuery(0),
    _needNewRHSJacobian(true),
    _needNewLHSJacobian(true)
{ // constructor
} // constructor

// ----------------------------------------------------------------------
// Destructor
pylith::feassemble::IntegratorPointwise::~IntegratorPointwise(void)
{ // destructor
    deallocate();
} // destructor

// ----------------------------------------------------------------------
// Deallocate PETSc and local data structures.
void
pylith::feassemble::IntegratorPointwise::deallocate(void)
{ // deallocate
    PYLITH_METHOD_BEGIN;

    delete _normalizer; _normalizer = 0;
    delete _logger; _logger = 0;
    delete _auxFields; _auxFields = 0;
    delete _auxFieldsQuery; _auxFieldsQuery = NULL;

    _auxFieldsDB = 0; // :TODO: Use shared pointer.

    PYLITH_METHOD_END;
} // deallocate

// ----------------------------------------------------------------------
// Return auxiliary fields for this problem
const pylith::topology::Field&
pylith::feassemble::IntegratorPointwise::auxFields(void) const
{ // auxFields
    PYLITH_METHOD_BEGIN;

    assert(_auxFields);

    PYLITH_METHOD_RETURN(*_auxFields);
} // auxFields

// ----------------------------------------------------------------------
// Check whether material has a given auxilirary field.
bool
pylith::feassemble::IntegratorPointwise::hasAuxField(const char* name)
{ // hasAuxField
    PYLITH_METHOD_BEGIN;

    assert(_auxFields);

    PYLITH_METHOD_RETURN(_auxFields->hasSubfield(name));
} // hasAuxField


// ----------------------------------------------------------------------
// Get auxiliary field.
void
pylith::feassemble::IntegratorPointwise::getAuxField(pylith::topology::Field *field,
                                                     const char* name) const
{ // getAuxField
    PYLITH_METHOD_BEGIN;

    journal::debug_t debug("integrator");
    debug << journal::at(__HERE__)
          << "IntegratorPointwise::getAuxField(field="<<field<<", name="<<name<<")" << journal::endl;

    assert(field);
    assert(_auxFields);

    field->copySubfield(*_auxFields, name);

    PYLITH_METHOD_END;
} // getAuxField


// ----------------------------------------------------------------------
// Set database for auxiliary fields.
void
pylith::feassemble::IntegratorPointwise::auxFieldsDB(spatialdata::spatialdb::SpatialDB* value) {
    _auxFieldsDB = value;
}

// ----------------------------------------------------------------------
// Set discretization information for auxiliary subfield.
void
pylith::feassemble::IntegratorPointwise::auxFieldDiscretization(const char* name,
                                                                const pylith::topology::FieldBase::DiscretizeInfo& feInfo)
{ // discretization
    journal::debug_t debug("integrator");
    debug << journal::at(__HERE__)
          << "IntegratorPointwise::auxFieldDiscretization(name="<<name<<", feInfo="<<&feInfo<<")" << journal::endl;

    _auxFieldsFEInfo[name] = feInfo;
} // discretization


// ----------------------------------------------------------------------
// Check whether RHS Jacobian needs to be recomputed.
bool
pylith::feassemble::IntegratorPointwise::needNewRHSJacobian(void) const {
    return _needNewRHSJacobian;
} // needNewRHSJacobian

// ----------------------------------------------------------------------
// Check whether LHS Jacobian needs to be recomputed.
bool
pylith::feassemble::IntegratorPointwise::needNewLHSJacobian(void) const {
    return _needNewLHSJacobian;
} // needNewLHSJacobian

// ----------------------------------------------------------------------
// Get discretization information for auxiliary subfield.
const pylith::topology::FieldBase::DiscretizeInfo&
pylith::feassemble::IntegratorPointwise::auxFieldDiscretization(const char* name) const
{ // discretization
    PYLITH_METHOD_BEGIN;

    discretizations_type::const_iterator iter = _auxFieldsFEInfo.find(name);
    if (iter != _auxFieldsFEInfo.end()) {
        PYLITH_METHOD_RETURN(iter->second);
    } else { // not found so try default
        iter = _auxFieldsFEInfo.find("default");
        if (iter == _auxFieldsFEInfo.end()) {
            throw std::logic_error("Default discretization not set for auxiliary fields.");
        } // if
    } // if/else

    PYLITH_METHOD_RETURN(iter->second); // default
} // discretization


// ----------------------------------------------------------------------
// Verify configuration is acceptable.
void
pylith::feassemble::IntegratorPointwise::verifyConfiguration(const pylith::topology::Field& solution) const
{ // verifyConfiguration
    PYLITH_METHOD_BEGIN;

    journal::debug_t debug("integrator");
    debug << journal::at(__HERE__)
          << "IntegratorPointwise::verifyConfiguration(solution="<<solution.label()<<") empty implementation." << journal::endl;

    PYLITH_METHOD_END;
} // verifyConfiguration


// ----------------------------------------------------------------------
// Set manager of scales used to nondimensionalize problem.
void
pylith::feassemble::IntegratorPointwise::normalizer(const spatialdata::units::Nondimensional& dim)
{ // normalizer
    journal::debug_t debug("integrator");
    debug << journal::at(__HERE__)
          << "IntegratorPointwise::normalizer(dim="<<&dim<<")" << journal::endl;

    if (0 == _normalizer)
        _normalizer = new spatialdata::units::Nondimensional(dim);
    else
        *_normalizer = dim;
} // normalizer


// ----------------------------------------------------------------------
// Update state variables as needed.
void
pylith::feassemble::IntegratorPointwise::updateStateVars(const pylith::topology::Field& solution)
{ // updateState
    PYLITH_METHOD_BEGIN;

    journal::debug_t debug("integrator");
    debug << journal::at(__HERE__)
          << "IntegratorPointwise::updateStateVars(solution) empty implementation." << journal::endl;

    PYLITH_METHOD_END;
} // updateStateVars



// End of file
