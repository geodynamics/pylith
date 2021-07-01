// -*- C++ -*-
//
// ----------------------------------------------------------------------
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
// ----------------------------------------------------------------------
//

#include <portinfo>

#include "DirichletUserFn.hh" // implementation of object methods

#include "pylith/feassemble/ConstraintUserFn.hh" // USES ConstraintBoundary
#include "pylith/topology/Field.hh" // USES Field
#include "pylith/topology/FieldOps.hh" // USES FieldOps
#include "pylith/topology/Mesh.hh" // USES Mesh

#include "spatialdata/units/Nondimensional.hh" // USES Nondimensional

#include "pylith/utils/error.hh" // USES PYLITH_METHOD_BEGIN/END
#include "pylith/utils/journals.hh" // USES PYLITH_COMPONENT_*

#include <cassert> // USES assert()
#include <stdexcept> // USES std::runtime_error
#include <sstream> // USES std::ostringstream
#include <typeinfo> // USES typeid()

// ---------------------------------------------------------------------------------------------------------------------
// Default constructor.
pylith::bc::DirichletUserFn::DirichletUserFn(void) :
    _fn(NULL),
    _fnDot(NULL) {
    PyreComponent::setName("dirichletuserfn");
} // constructor


// ---------------------------------------------------------------------------------------------------------------------
// Destructor.
pylith::bc::DirichletUserFn::~DirichletUserFn(void) {
    deallocate();
} // destructor


// ---------------------------------------------------------------------------------------------------------------------
// Deallocate PETSc and local data structures.
void
pylith::bc::DirichletUserFn::deallocate(void) {
    PYLITH_METHOD_BEGIN;

    BoundaryCondition::deallocate();

    PYLITH_METHOD_END;
} // deallocate


// ---------------------------------------------------------------------------------------------------------------------
// Set indices of constrained degrees of freedom at each location.
void
pylith::bc::DirichletUserFn::setConstrainedDOF(const int* flags,
                                               const int size) {
    PYLITH_COMPONENT_DEBUG("setConstrainedDOF(flags="<<flags<<", size"<<size<<")");

    assert((flags && size > 0) || (!flags && 0 == size) );
    _constrainedDOF.resize(size);
    for (int i = 0; i < size; ++i) {
        _constrainedDOF[i] = flags[i];
    } // for
} // setConstrainedDOF


// ---------------------------------------------------------------------------------------------------------------------
// Get indices of constrained degrees of freedom.
const pylith::int_array&
pylith::bc::DirichletUserFn::getConstrainedDOF(void) const {
    return _constrainedDOF;
} // getConstrainedDOF


// ---------------------------------------------------------------------------------------------------------------------
// Set user function specifying field on boundary.
void
pylith::bc::DirichletUserFn::setUserFn(PetscUserFieldFunc fn) {
    PYLITH_COMPONENT_DEBUG("setUserFn(fn="<<fn<<")");

    _fn = fn;
} // setUserFn


// ---------------------------------------------------------------------------------------------------------------------
// Get user function specifying field on boundary.
PetscUserFieldFunc
pylith::bc::DirichletUserFn::getUserFn(void) const {
    return _fn;
} // getUserFn


// ---------------------------------------------------------------------------------------------------------------------
// Set user function specifying time derivative of field on boundary.
void
pylith::bc::DirichletUserFn::setUserFnDot(PetscUserFieldFunc fn) {
    PYLITH_COMPONENT_DEBUG("setUserFnDot(fn="<<fn<<")");

    _fnDot = fn;
} // setUserFnDot


// ---------------------------------------------------------------------------------------------------------------------
// Get user function specifying time derivative of field on boundary.
PetscUserFieldFunc
pylith::bc::DirichletUserFn::getUserFnDot(void) const {
    return _fnDot;
} // getUserFnDot


// ---------------------------------------------------------------------------------------------------------------------
// Verify configuration is acceptable.
void
pylith::bc::DirichletUserFn::verifyConfiguration(const pylith::topology::Field& solution) const {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("verifyConfiguration(solution="<<solution.getLabel()<<")");

    if (!solution.hasSubfield(_subfieldName.c_str())) {
        std::ostringstream msg;
        msg << "Cannot constrain field '"<< _subfieldName
            << "' in component '" << PyreComponent::getIdentifier() << "'"
            << "; field is not in solution.";
        throw std::runtime_error(msg.str());
    } // if

    const topology::Field::SubfieldInfo& info = solution.getSubfieldInfo(_subfieldName.c_str());
    const int numComponents = info.description.numComponents;
    const int numConstrained = _constrainedDOF.size();
    for (int iConstrained = 0; iConstrained < numConstrained; ++iConstrained) {
        if (_constrainedDOF[iConstrained] >= numComponents) {
            std::ostringstream msg;
            msg << "Cannot constrain degree of freedom '" << _constrainedDOF[iConstrained] << "'"
                << " in component '" << PyreComponent::getIdentifier() << "'"
                << "; solution field '" << _subfieldName << "' contains only " << numComponents << " components.";
            throw std::runtime_error(msg.str());
        } // if
    } // for

    PYLITH_METHOD_END;
} // verifyConfiguration


// ---------------------------------------------------------------------------------------------------------------------
// Create integrator and set kernels.
pylith::feassemble::Integrator*
pylith::bc::DirichletUserFn::createIntegrator(const pylith::topology::Field& solution) {
    PYLITH_COMPONENT_DEBUG("createIntegrator(solution="<<solution.getLabel()<<") empty method");

    return NULL;
} // createIntegrator


// ---------------------------------------------------------------------------------------------------------------------
// Create constraint and set kernels.
pylith::feassemble::Constraint*
pylith::bc::DirichletUserFn::createConstraint(const pylith::topology::Field& solution) {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("createConstraint(solution="<<solution.getLabel()<<")");

    pylith::feassemble::ConstraintUserFn* constraint = new pylith::feassemble::ConstraintUserFn(this);assert(constraint);
    constraint->setMarkerLabel(getMarkerLabel());
    constraint->setConstrainedDOF(&_constrainedDOF[0], _constrainedDOF.size());
    constraint->setSubfieldName(_subfieldName.c_str());
    constraint->setUserFn(_fn);
    constraint->setUserFnDot(_fnDot);

    PYLITH_METHOD_RETURN(constraint);
} // createConstraint


// ---------------------------------------------------------------------------------------------------------------------
// Create auxiliary field.
pylith::topology::Field*
pylith::bc::DirichletUserFn::createAuxiliaryField(const pylith::topology::Field& solution,
                                                  const pylith::topology::Mesh& domainMesh) {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("createAuxiliaryField(solution="<<solution.getLabel()<<", domainMesh=)"<<typeid(domainMesh).name()<<") empty method");

    PYLITH_METHOD_RETURN(NULL);
} // createAuxiliaryField


// ---------------------------------------------------------------------------------------------------------------------
// Create derived field.
pylith::topology::Field*
pylith::bc::DirichletUserFn::createDerivedField(const pylith::topology::Field& solution,
                                                const pylith::topology::Mesh& domainMesh) {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("createDerivedField(solution="<<solution.getLabel()<<", domainMesh=)"<<typeid(domainMesh).name()<<") empty method");

    PYLITH_METHOD_RETURN(NULL);
} // createDerivedField


// ---------------------------------------------------------------------------------------------------------------------
// Get auxiliary factory associated with physics.
pylith::feassemble::AuxiliaryFactory*
pylith::bc::DirichletUserFn::_getAuxiliaryFactory(void) {
    return NULL;
} // _getAuxiliaryFactory


// End of file
