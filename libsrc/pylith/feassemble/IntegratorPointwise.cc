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

#include "pylith/feassemble/AuxiliaryFactory.hh" // USES AuxiliaryFactory
#include "pylith/topology/Mesh.hh" // USES Mesh
#include "pylith/topology/Field.hh" // USES Field
#include "pylith/meshio/OutputManager.hh" // USES OutputManager

#include "spatialdata/units/Nondimensional.hh" // USES Nondimensional

#include "pylith/utils/EventLogger.hh" // USES EventLogger
#include "pylith/utils/journals.hh" // USES PYLITH_COMPONENT_*

#include <cassert> // USES assert()
#include <typeinfo> // USES typeid()
#include <stdexcept> \
    // USES std::runtime_error

// ----------------------------------------------------------------------
// Constructor
pylith::feassemble::IntegratorPointwise::IntegratorPointwise(void) :
    _normalizer(new spatialdata::units::Nondimensional),
    _gravityField(NULL),
    _auxField(NULL),
    _output(NULL),
    _logger(NULL),
    _needNewRHSJacobian(true),
    _needNewLHSJacobian(true)
{} // constructor

// ----------------------------------------------------------------------
// Destructor
pylith::feassemble::IntegratorPointwise::~IntegratorPointwise(void) {
    deallocate();
} // destructor

// ----------------------------------------------------------------------
// Deallocate PETSc and local data structures.
void
pylith::feassemble::IntegratorPointwise::deallocate(void) {
    PYLITH_METHOD_BEGIN;

    delete _normalizer; _normalizer = NULL;
    delete _logger; _logger = NULL;
    delete _auxField; _auxField = NULL;

    _gravityField = NULL; // :KLUDGE: Memory managed by Python object. :TODO: Use shared pointer.
    _output = NULL; // :KLUDGE: Memory managed by Python object. :TODO: Use shared pointer.

    PYLITH_METHOD_END;
} // deallocate

// ----------------------------------------------------------------------
// Get auxiliary field.
const pylith::topology::Field&
pylith::feassemble::IntegratorPointwise::auxField(void) const {
    PYLITH_METHOD_BEGIN;

    assert(_auxField);

    PYLITH_METHOD_RETURN(*_auxField);
} // auxField

// ----------------------------------------------------------------------
// Set database for auxiliary fields.
void
pylith::feassemble::IntegratorPointwise::auxFieldDB(spatialdata::spatialdb::SpatialDB* value) {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("auxFieldDB(value="<<value<<")");

    pylith::feassemble::AuxiliaryFactory* factory = _auxFactory(); assert(factory);
    factory->queryDB(value);

    PYLITH_METHOD_END;
} // auxFieldDB


// ----------------------------------------------------------------------
// Set discretization information for auxiliary subfield.
void
pylith::feassemble::IntegratorPointwise::auxSubfieldDiscretization(const char* name,
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
pylith::feassemble::IntegratorPointwise::normalizer(const spatialdata::units::Nondimensional& dim) {
    PYLITH_COMPONENT_DEBUG("normalizer(dim="<<typeid(dim).name()<<")");

    if (!_normalizer) {
        _normalizer = new spatialdata::units::Nondimensional(dim);
    } else {
        *_normalizer = dim;
    } // if/else
} // normalizer


// ----------------------------------------------------------------------
// Set gravity field.
void
pylith::feassemble::IntegratorPointwise::gravityField(spatialdata::spatialdb::GravityField* const g) {
    _gravityField = g;
} // gravityField


// ----------------------------------------------------------------------
// Set output manager.
void
pylith::feassemble::IntegratorPointwise::output(pylith::meshio::OutputManager* manager) {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("output(manager="<<manager<<")");

    _output = manager;

    PYLITH_METHOD_END;
} // output


// ----------------------------------------------------------------------
// Update auxiliary fields at beginning of time step.
void
pylith::feassemble::IntegratorPointwise::prestep(const double t,
                                                 const double dt) {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("prestep(t="<<t<<", dt="<<dt<<") empty method");

    // Default is to do nothing.

    PYLITH_METHOD_END;
} // prestep


// ----------------------------------------------------------------------
// Update state variables as needed.
void
pylith::feassemble::IntegratorPointwise::updateStateVars(const pylith::topology::Field& solution) {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("updateStateVars(solution="<<solution.label()<<")");

    PYLITH_COMPONENT_ERROR(":TODO: @brad Implement updateStateVars().");

    PYLITH_METHOD_END;
} // updateStateVars


// ----------------------------------------------------------------------
// Write information (auxiliary field) output.
void
pylith::feassemble::IntegratorPointwise::writeInfo(void) {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("writeInfo(void)");

    if (_output) {
        assert(_auxField);
        _output->writeInfo(*_auxField);
    } // if

    PYLITH_METHOD_END;
} // writeInfo


// ----------------------------------------------------------------------
// Write solution related output.
void
pylith::feassemble::IntegratorPointwise::writeTimeStep(const PylithReal t,
                                                       const PylithInt tindex,
                                                       const pylith::topology::Field& solution) {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("writeTimeStep(t="<<t<<", tindex="<<tindex<<", solution="<<solution.label()<<")");

    if (_output) {
        assert(_auxField);
        _output->writeTimeStep(t, tindex, solution, *_auxField);
    } // if

    PYLITH_METHOD_END;
} // writeTimeStep


// End of file
