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
// Copyright (c) 2010-2016 University of California, Davis
//
// See COPYING for license information.
//
// ----------------------------------------------------------------------
//

#include <portinfo>

#include "ConstraintPointwise.hh" // implementation of object methods

#include "spatialdata/units/Nondimensional.hh" // USES Nondimensional

#include "pylith/topology/Mesh.hh" // USES Mesh
#include "pylith/topology/Field.hh" // USES Field
#include "pylith/topology/FieldQuery.hh" // USES FieldQuery

#include "pylith/utils/journals.hh" // USES PYLITH_JOURNAL_*

#include "spatialdata/units/Nondimensional.hh" // USES Nondimensional

#include "pylith/utils/EventLogger.hh" // USES EventLogger


#include <cassert> // USES assert()
#include <stdexcept> // USES std::runtime_error

// ----------------------------------------------------------------------
// Default constructor.
pylith::feassemble::ConstraintPointwise::ConstraintPointwise(void) :
    _normalizer(new spatialdata::units::Nondimensional),
    _logger(0),
    _auxFields(0),
    _auxFieldsDB(0),
    _auxFieldsQuery(0)
{ // constructor
} // constructor

// ----------------------------------------------------------------------
// Destructor.
pylith::feassemble::ConstraintPointwise::~ConstraintPointwise(void)
{ // destructor
    deallocate();
} // destructor

// ----------------------------------------------------------------------
// Deallocate PETSc and local data structures.
void
pylith::feassemble::ConstraintPointwise::deallocate(void)
{ // deallocate
    PYLITH_METHOD_BEGIN;

    delete _normalizer; _normalizer = 0;
    delete _logger; _logger = 0;
    delete _auxFields; _auxFields = 0;
    delete _auxFieldsQuery; _auxFieldsQuery = 0;

    _auxFieldsDB = 0; // :TODO: Use shared pointer.

    PYLITH_METHOD_END;
} // deallocate


// ----------------------------------------------------------------------
// Set name of field in solution to constrain.
void
pylith::feassemble::ConstraintPointwise::field(const char* value)
{  // field
    PYLITH_METHOD_BEGIN;

    if (strlen(value) == 0) {
        throw std::runtime_error("Empty string given for name of solution field to constrain.");
    } // if
    _field = value;

    PYLITH_METHOD_END;
}  // field


// ----------------------------------------------------------------------
// Get name of field in solution to constrain.
const char*
pylith::feassemble::ConstraintPointwise::field(void) const
{ // field
    journal::debug_t debug("constraint");
    debug << journal::at(__HERE__)
          << "ConstraintPointwise::field()" << journal::endl;

    return _field.c_str();
} // field


// ----------------------------------------------------------------------
// Set indices of constrained degrees of freedom at each location.
void
pylith::feassemble::ConstraintPointwise::constrainedDOF(const int* flags,
                                                        const int size)
{ // constrainedDOF
    PYLITH_METHOD_BEGIN;
    PYLITH_JOURNAL_DEBUG("constrainedDOF(flags="<<flags<<", size="<<size<<")");

    assert((size > 0 && flags) || (!size && !flags));

    _constrainedDOF.resize(size);
    for (int i=0; i < size; ++i) {
        if (flags[i] < 0 || flags[i] > 2) {
            std::ostringstream msg;
            msg << "Constrained DOF '" << flags[i] << "' out of range"
                << " in component '" << PyreComponent::identifier() << "'"
                << ". DOF must be in range [0,2].";
            throw std::runtime_error(msg.str());
        } // if
        _constrainedDOF[i] = flags[i];
    } // for

    PYLITH_METHOD_END;
} // constrainedDOF


// ----------------------------------------------------------------------
// Get indices of constrained degrees of freedom.
const pylith::int_array&
pylith::feassemble::ConstraintPointwise::constrainedDOF(void) const
{ // constrainedDOF
    return _constrainedDOF;
} // constrainedDOF

// ----------------------------------------------------------------------
// Return auxiliary fields for this problem
const pylith::topology::Field&
pylith::feassemble::ConstraintPointwise::auxFields(void) const
{ // auxFields
    PYLITH_METHOD_BEGIN;

    assert(_auxFields);

    PYLITH_METHOD_RETURN(*_auxFields);
} // auxFields

// ----------------------------------------------------------------------
// Check whether constraint has a given auxilirary field.
bool
pylith::feassemble::ConstraintPointwise::hasAuxField(const char* name)
{ // hasAuxField
    PYLITH_METHOD_BEGIN;

    assert(_auxFields);

    PYLITH_METHOD_RETURN(_auxFields->hasSubfield(name));
} // hasAuxField


// ----------------------------------------------------------------------
// Get auxiliary field.
void
pylith::feassemble::ConstraintPointwise::getAuxField(pylith::topology::Field *field,
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
pylith::feassemble::ConstraintPointwise::auxFieldsDB(spatialdata::spatialdb::SpatialDB* value) {
    _auxFieldsDB = value;
}

// ----------------------------------------------------------------------
// Set discretization information for auxiliary subfield.
void
pylith::feassemble::ConstraintPointwise::auxFieldDiscretization(const char* name,
                                                                const int basisOrder,
                                                                const int quadOrder,
                                                                const bool isBasisContinuous)
{ // discretization
    PYLITH_JOURNAL_DEBUG("auxFieldDiscretization(name="<<name<<", basisOrder="<<basisOrder<<", quadOrder="<<quadOrder<<", isBasisContinuous="<<isBasisContinuous<<")");

    pylith::topology::FieldBase::DiscretizeInfo feInfo;
    feInfo.basisOrder = basisOrder;
    feInfo.quadOrder = quadOrder;
    feInfo.isBasisContinuous = isBasisContinuous;
    _auxFieldsFEInfo[name] = feInfo;
} // discretization


// ----------------------------------------------------------------------
// Get discretization information for auxiliary subfield.
const pylith::topology::FieldBase::DiscretizeInfo&
pylith::feassemble::ConstraintPointwise::auxFieldDiscretization(const char* name) const
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
// Set manager of scales used to nondimensionalize problem.
void
pylith::feassemble::ConstraintPointwise::normalizer(const spatialdata::units::Nondimensional& dim)
{ // normalizer
    PYLITH_JOURNAL_DEBUG("normalizer(dim="<<&dim<<")");

    if (!_normalizer) {
        _normalizer = new spatialdata::units::Nondimensional(dim);
    } else {
        *_normalizer = dim;
    } // if/else
} // normalizer


// ----------------------------------------------------------------------
// Verify configuration is acceptable.
void
pylith::feassemble::ConstraintPointwise::verifyConfiguration(const pylith::topology::Field& solution) const
{ // verifyConfiguration
    PYLITH_METHOD_BEGIN;
    PYLITH_JOURNAL_DEBUG("verifyConfiguration(solution="<<solution.label()<<")");

    if (!solution.hasSubfield(_field.c_str())) {
        std::ostringstream msg;
        msg << "Cannot constrain field '"<< _field
            << "' in component '" << PyreComponent::identifier() << "'"
            << "; field is not in solution.";
        throw std::runtime_error(msg.str());
    } // if

    const topology::Field::SubfieldInfo& info = solution.subfieldInfo(_field.c_str());
    const int numComponents = info.numComponents;
    const int numConstrained = _constrainedDOF.size();
    for (int iConstrained=0; iConstrained < numConstrained; ++iConstrained) {
        if (_constrainedDOF[iConstrained] >= numComponents) {
            std::ostringstream msg;
            msg << "Cannot constrain degree of freedom '" << _constrainedDOF[iConstrained] << "'"
                << " in component '" << PyreComponent::identifier() << "'"
                << "; solution field '" << _field << "' contains only " << numComponents << " components.";
            throw std::runtime_error(msg.str());
        } // if
    } // for

    PYLITH_METHOD_END;
} // verifyConfiguration


// ----------------------------------------------------------------------
// Update auxiliary fields at beginning of time step.
void
pylith::feassemble::ConstraintPointwise::prestep(const double t)
{ // prestep
  // Default is to do nothing.
} // prestep


// End of file
