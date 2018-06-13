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

#include "pylith/feassemble/ConstraintPointwise.hh" // implementation of object methods

#include "pylith/feassemble/AuxiliaryFactory.hh" // USES AuxiliaryFactory
#include "pylith/topology/Mesh.hh" // USES Mesh
#include "pylith/topology/Field.hh" // USES Field

#include "spatialdata/units/Nondimensional.hh" // USES Nondimensional
#include "spatialdata/geocoords/CoordSys.hh" // USES CoordSys

#include "pylith/utils/journals.hh" // USES PYLITH_COMPONENT_*
#include "pylith/utils/EventLogger.hh" // USES EventLogger

#include <cassert> // USES assert()
#include <typeinfo> // USES typeid()
#include <stdexcept> \
    // USES std::runtime_error

// ----------------------------------------------------------------------
// Default constructor.
pylith::feassemble::ConstraintPointwise::ConstraintPointwise(void) :
    _normalizer(new spatialdata::units::Nondimensional),
    _auxField(NULL),
    _logger(NULL)
{}

// ----------------------------------------------------------------------
// Destructor.
pylith::feassemble::ConstraintPointwise::~ConstraintPointwise(void) {
    deallocate();
} // destructor

// ----------------------------------------------------------------------
// Deallocate PETSc and local data structures.
void
pylith::feassemble::ConstraintPointwise::deallocate(void) {
    PYLITH_METHOD_BEGIN;

    ObservedComponent::deallocate();

    delete _normalizer; _normalizer = NULL;
    delete _logger; _logger = NULL;
    delete _auxField; _auxField = NULL;

    PYLITH_METHOD_END;
} // deallocate


// ----------------------------------------------------------------------
// Set indices of constrained degrees of freedom at each location.
void
pylith::feassemble::ConstraintPointwise::constrainedDOF(const int* flags,
                                                        const int size) {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("constrainedDOF(flags="<<flags<<", size="<<size<<")");

    assert((size > 0 && flags) || (!size && !flags));

    _constrainedDOF.resize(size);
    for (int i = 0; i < size; ++i) {
        if (flags[i] < 0) {
            std::ostringstream msg;
            msg << "Constrained DOF '" << flags[i] << "' must be nonnegative in component '" << PyreComponent::identifier() << "'.";
            throw std::runtime_error(msg.str());
        } // if
        _constrainedDOF[i] = flags[i];
    } // for

    PYLITH_METHOD_END;
} // constrainedDOF


// ----------------------------------------------------------------------
// Get indices of constrained degrees of freedom.
const pylith::int_array&
pylith::feassemble::ConstraintPointwise::constrainedDOF(void) const {
    return _constrainedDOF;
} // constrainedDOF

// ----------------------------------------------------------------------
// Return auxiliary subfields for this problem.
const pylith::topology::Field&
pylith::feassemble::ConstraintPointwise::auxField(void) const {
    PYLITH_METHOD_BEGIN;

    assert(_auxField);
    PYLITH_METHOD_RETURN(*_auxField);
} // auxField

// ----------------------------------------------------------------------
// Set database for filling auxiliary subfields.
void
pylith::feassemble::ConstraintPointwise::auxFieldDB(spatialdata::spatialdb::SpatialDB* value) {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("auxFieldDB(value="<<value<<")");

    pylith::feassemble::AuxiliaryFactory* factory = _auxFactory(); assert(factory);
    factory->queryDB(value);

    PYLITH_METHOD_END;
} // auxFieldDB

// ----------------------------------------------------------------------
// Set discretization information for auxiliary subfield.
void
pylith::feassemble::ConstraintPointwise::auxSubfieldDiscretization(const char* name,
                                                                   const int basisOrder,
                                                                   const int quadOrder,
                                                                   const bool isBasisContinuous,
                                                                   const pylith::topology::FieldBase::SpaceEnum feSpace) {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("auxSubfieldDiscretization(name="<<name<<", basisOrder="<<basisOrder<<", quadOrder="<<quadOrder<<", isBasisContinuous="<<isBasisContinuous<<")");

    pylith::feassemble::AuxiliaryFactory* factory = _auxFactory(); assert(factory);
    factory->subfieldDiscretization(name, basisOrder, quadOrder, isBasisContinuous, feSpace);

    PYLITH_METHOD_END;
} // auxSubfieldDiscretization


// ----------------------------------------------------------------------
// Set manager of scales used to nondimensionalize problem.
void
pylith::feassemble::ConstraintPointwise::normalizer(const spatialdata::units::Nondimensional& dim) {
    PYLITH_COMPONENT_DEBUG("normalizer(dim="<<typeid(dim).name()<<")");

    if (!_normalizer) {
        _normalizer = new spatialdata::units::Nondimensional(dim);
    } else {
        *_normalizer = dim;
    } // if/else
} // normalizer


// ----------------------------------------------------------------------
// Update auxiliary fields at beginning of time step.
void
pylith::feassemble::ConstraintPointwise::prestep(const double t,
                                                 const double dt)
{ // prestep
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("prestep(t="<<t<<", dt="<<dt<<") empty method");

    // Default is to do nothing.

    PYLITH_METHOD_END;
} // prestep


// End of file
