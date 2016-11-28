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

#include "spatialdata/units/Nondimensional.hh" // USES Nondimensional

#include "pylith/utils/EventLogger.hh" // USES EventLogger

#include "journal/debug.h" // USES journal::debug_t

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

    journal::debug_t debug("constraint");
    debug << journal::at(__HERE__)
          << "ConstraintPointwise::constrainedDOF(flags="<<flags<<", size="<<size<<")" << journal::endl;

    if (size > 0) {
        assert(flags);
    } // if

    _constrainedDOF.resize(size);
    for (int i=0; i < size; ++i) {
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

    journal::debug_t debug("constraint");
    debug << journal::at(__HERE__)
          << "ConstraintPointwise::getAuxField(field="<<field<<", name="<<name<<")" << journal::endl;

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
    journal::debug_t debug("constraint");
    debug << journal::at(__HERE__)
          << "ConstraintPointwise::auxFieldDiscretization(name="<<name<<", basisOrder="<<basisOrder<<", quadOrder="<<quadOrder<<", isBasisContinuous="<<isBasisContinuous<<")" << journal::endl;

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
    journal::debug_t debug("constraint");
    debug << journal::at(__HERE__)
          << "ConstraintPointwise::normalizer(dim="<<&dim<<")" << journal::endl;

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

    journal::debug_t debug("constraint");
    debug << journal::at(__HERE__)
          << "ConstraintPointwise::verifyConfiguration(solution="<<solution.label()<<") empty implementation." << journal::endl;

    PYLITH_METHOD_END;
} // verifyConfiguration


// ----------------------------------------------------------------------
// Update auxiliary fields at beginning of time step.
void
pylith::feassemble::ConstraintPointwise::prestep(const double t)
{ // prestep

} // prestep


// End of file
