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
#include "pylith/utils/journals.hh" // USES PYLITH_JOURNAL_*

#include <cassert> // USES assert()
#include <stdexcept> // USES std::runtime_error

// ----------------------------------------------------------------------
// Constructor
pylith::feassemble::IntegratorPointwise::IntegratorPointwise(void) :
    _normalizer(new spatialdata::units::Nondimensional),
    _gravityField(0),
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
    delete _auxFieldsQuery; _auxFieldsQuery = 0;

    _gravityField = 0; // :TODO: Use shared points.
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
    PYLITH_JOURNAL_DEBUG("getAuxField(field="<<field<<", name="<<name<<")");

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
                                                                const int basisOrder,
                                                                const int quadOrder,
                                                                const bool isBasisContinuous)
{ // discretization
    PYLITH_METHOD_BEGIN;
    PYLITH_JOURNAL_DEBUG("auxFieldDiscretization(name="<<name<<", basisOrder="<<basisOrder<<", quadOrder="<<quadOrder<<", isBasisContinuous="<<isBasisContinuous<<")");

    pylith::topology::FieldBase::DiscretizeInfo feInfo;
    feInfo.basisOrder = basisOrder;
    feInfo.quadOrder = quadOrder;
    feInfo.isBasisContinuous = isBasisContinuous;
    _auxFieldsFEInfo[name] = feInfo;

    PYLITH_METHOD_END;
} // discretization


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
// Set manager of scales used to nondimensionalize problem.
void
pylith::feassemble::IntegratorPointwise::normalizer(const spatialdata::units::Nondimensional& dim)
{ // normalizer
    PYLITH_JOURNAL_DEBUG("normalizer(dim="<<&dim<<")");

    if (0 == _normalizer) {
        _normalizer = new spatialdata::units::Nondimensional(dim);
    } else {
        *_normalizer = dim;
    } // if/else
} // normalizer


// ----------------------------------------------------------------------
// Set gravity field.
void
pylith::feassemble::IntegratorPointwise::gravityField(spatialdata::spatialdb::GravityField* const g)
{ // gravityField
    _gravityField = g;
} // gravityField


// ----------------------------------------------------------------------
// Update state variables as needed.
void
pylith::feassemble::IntegratorPointwise::updateStateVars(const pylith::topology::Field& solution)
{ // updateState
    PYLITH_METHOD_BEGIN;
    PYLITH_JOURNAL_DEBUG("updateStateVars(solution="<<solution.label()<<")");

    PYLITH_JOURNAL_ERROR(":TODO: @brad Implement updateStateVars().");

    PYLITH_METHOD_END;
} // updateStateVars


// ----------------------------------------------------------------------
// Update auxiliary fields at beginning of time step.
void
pylith::feassemble::IntegratorPointwise::prestep(const double t,
                                                 const double dt)
{ // prestep
  // Default is to do nothing.
} // prestep


// End of file
