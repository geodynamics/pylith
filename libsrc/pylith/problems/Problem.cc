// -*- C++ -*-
//
// ======================================================================
//
// Brad T. Aagaard, U.S. Geological Survey
// Charles A. Williams, GNS Science
// Matthew G. Knepley, University of Chicago
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

#include "Problem.hh" // implementation of class methods

#include "pylith/topology/Mesh.hh" // USES Mesh
#include "pylith/topology/Field.hh" // USES Field

#include "pylith/feassemble/IntegratorPointwise.hh" // USES IntegratorPointwise
#include "pylith/feassemble/ConstraintPointwise.hh" // USES ConstraintPointwise
#include "pylith/topology/Jacobian.hh" // USES Jacobian
#include "pylith/topology/MeshOps.hh" // USES MeshOps

#include "spatialdata/units/Nondimensional.hh" // USES Nondimensional

#include "pylith/utils/error.hh" // USES PYLITH_CHECK_ERROR
#include "pylith/utils/journals.hh" // USES PYLITH_JOURNAL_*

#include <cassert> // USES assert()

// ----------------------------------------------------------------------
// Constructor
pylith::problems::Problem::Problem() :
    _solution(0),
    _jacobianLHSLumpedInv(0),
    _normalizer(0),
    _integrators(0),
    _constraints(0),
    _outputs(0),
    _solverType(LINEAR)
{ // constructor
} // constructor

// ----------------------------------------------------------------------
// Destructor
pylith::problems::Problem::~Problem(void)
{ // destructor
    deallocate();
} // destructor

// ----------------------------------------------------------------------
// Deallocate PETSc and local data structures.
void
pylith::problems::Problem::deallocate(void)
{ // deallocate
    PYLITH_METHOD_BEGIN;

    _solution = 0; // Held by Python. :KLUDGE: :TODO: Use shared pointer.
    delete _jacobianLHSLumpedInv; _jacobianLHSLumpedInv = 0;
    delete _normalizer; _normalizer = 0;

    PYLITH_METHOD_END;
} // deallocate

// ----------------------------------------------------------------------
// Set problem type.
void
pylith::problems::Problem::solverType(const SolverTypeEnum value)
{ // solverType
    PYLITH_JOURNAL_DEBUG("Problem::solverType(value="<<value<<")");

    _solverType = value;
} // solverType

// ----------------------------------------------------------------------
// Get problem type.
pylith::problems::Problem::SolverTypeEnum
pylith::problems::Problem::solverType(void) const
{ // solverType
    return _solverType;
} // solverType

// ----------------------------------------------------------------------
// Set manager of scales used to nondimensionalize problem.
void
pylith::problems::Problem::normalizer(const spatialdata::units::Nondimensional& dim)
{ // normalizer
    PYLITH_JOURNAL_DEBUG("Problem::normalizer(dim="<<&dim<<")");

    if (!_normalizer) {
        _normalizer = new spatialdata::units::Nondimensional(dim);
    } else {
        *_normalizer = dim;
    } // if/else
} // normalizer

// ----------------------------------------------------------------------
// Set gravity field.
void
pylith::problems::Problem::gravityField(spatialdata::spatialdb::GravityField* const g)
{ // gravityField
    _gravityField = g;
} // gravityField

// ----------------------------------------------------------------------
// Set solution field.
void
pylith::problems::Problem::solution(pylith::topology::Field* field)
{ // solution
    _solution = field;
} // solution

// ----------------------------------------------------------------------
// Set integrators over the mesh.
void
pylith::problems::Problem::integrators(pylith::feassemble::IntegratorPointwise* integratorArray[],
                                       const int numIntegrators)
{ // integrators
    PYLITH_METHOD_BEGIN;
    PYLITH_JOURNAL_DEBUG("Problem::integrators("<<integratorArray<<", numIntegrators="<<numIntegrators<<")");

    assert( (!integratorArray && 0 == numIntegrators) || (integratorArray && 0 < numIntegrators) );

    _integrators.resize(numIntegrators);
    for (int i=0; i < numIntegrators; ++i) {
        _integrators[i] = integratorArray[i];
    } // for

    PYLITH_METHOD_END;
} // integrators

// ----------------------------------------------------------------------
// Set constraints over the mesh.
void
pylith::problems::Problem::constraints(pylith::feassemble::ConstraintPointwise* constraintArray[],
                                       const int numConstraints)
{ // constraints
    PYLITH_METHOD_BEGIN;
    PYLITH_JOURNAL_DEBUG("Problem::constraints("<<constraintArray<<", numConstraints="<<numConstraints<<")");

    assert( (!constraintArray && 0 == numConstraints) || (constraintArray && 0 < numConstraints) );

    _constraints.resize(numConstraints);
    for (int i=0; i < numConstraints; ++i) {
        _constraints[i] = constraintArray[i];
    } // for

    PYLITH_METHOD_END;
} // constraints

// ----------------------------------------------------------------------
// Set constraints over the mesh.
void
pylith::problems::Problem::outputs(pylith::meshio::OutputManager* outputArray[],
                                   const int numOutputs)
{ // outputs
    PYLITH_METHOD_BEGIN;
    PYLITH_JOURNAL_DEBUG("Problem::output("<<outputArray<<", numOutputs="<<numOutputs<<")");

    assert( (!outputArray && 0 == numOutputs) || (outputArray && 0 < numOutputs) );

    _outputs.resize(numOutputs);
    for (int i=0; i < numOutputs; ++i) {
        _outputs[i] = outputArray[i];
    } // for

    PYLITH_METHOD_END;
} // outputs


// ----------------------------------------------------------------------
// Do minimal initialization.
void
pylith::problems::Problem::preinitialize(const pylith::topology::Mesh& mesh)
{ // preinitialize
    PYLITH_METHOD_BEGIN;
    PYLITH_JOURNAL_DEBUG("Problem::preinitialzie(mesh="<<&mesh<<")");

    assert(_normalizer);

    const size_t numIntegrators = _integrators.size();
    for (size_t i=0; i < numIntegrators; ++i) {
        assert(_integrators[i]);
        _integrators[i]->normalizer(*_normalizer);
        _integrators[i]->gravityField(_gravityField);
    } // for

    const size_t numConstraints = _constraints.size();
    for (size_t i=0; i < numConstraints; ++i) {
        assert(_constraints[i]);
        _constraints[i]->normalizer(*_normalizer);
    } // for

    PYLITH_METHOD_END;
} // preinitialize

// ----------------------------------------------------------------------
// Verify configuration.
void
pylith::problems::Problem::verifyConfiguration(int* const materialIds,
                                               const int numMaterials)
{ // verifyConfiguration
    PYLITH_METHOD_BEGIN;
    PYLITH_JOURNAL_DEBUG("Problem::verifyConfiguration(materialIds="<<materialIds<<", numMaterials="<<numMaterials<<")");

    // Check to make sure material-id for each cell matches the id of a material.
    pylith::topology::MeshOps::checkMaterialIds(_solution->mesh(), materialIds, numMaterials);

    // Check to make sure integrators are compatible with the solution.
    const size_t numIntegrators = _integrators.size();
    for (size_t i=0; i < numIntegrators; ++i) {
        assert(_integrators[i]);
        _integrators[i]->verifyConfiguration(*_solution);
    } // for

    // Check to make sure constraints are compatible with the solution.
    const size_t numConstraints = _constraints.size();
    for (size_t i=0; i < numConstraints; ++i) {
        assert(_constraints[i]);
        _constraints[i]->verifyConfiguration(*_solution);
    } // for

    PYLITH_METHOD_END;
}  // verifyConfiguration

// ----------------------------------------------------------------------
// Initialize.
void
pylith::problems::Problem::initialize(void)
{ // initialize
    PYLITH_METHOD_BEGIN;
    PYLITH_JOURNAL_DEBUG("Problem::initialize()");

    assert(_solution);

    // Initialize integrators.
    const size_t numIntegrators = _integrators.size();
    for (size_t i=0; i < numIntegrators; ++i) {
        assert(_integrators[i]);
        _integrators[i]->initialize(*_solution);
    } // for

    // Initialize constraints.
    const size_t numConstraints = _constraints.size();
    for (size_t i=0; i < numConstraints; ++i) {
        assert(_constraints[i]);
        _constraints[i]->initialize(*_solution);
    } // for

    // Initialize solution field.
    _solution->allocate();
    _solution->zeroAll();

    PYLITH_METHOD_END;
} // initialize

// ----------------------------------------------------------------------
// Compute RHS residual for G(t,s).
void
pylith::problems::Problem::computeRHSResidual(PetscVec residualVec,
                                              const PylithReal t,
                                              const PylithReal dt,
                                              PetscVec solutionVec)
{ // computeRHSResidual
    PYLITH_METHOD_BEGIN;
    PYLITH_JOURNAL_DEBUG("Problem::computeRHSResidual(t="<<t<<", dt="<<dt<<", solutionVec="<<solutionVec<<", residualVec="<<residualVec<<")");

    assert(residualVec);
    assert(solutionVec);
    assert(_solution);

    // Update PyLith view of the solution.
    _solution->scatterGlobalToLocal(solutionVec);

    // Sum residual contributions across integrators.
    const size_t numIntegrators = _integrators.size();
    assert(numIntegrators > 0); // must have at least 1 integrator
    for (size_t i=0; i < numIntegrators; ++i) {
        _integrators[i]->computeRHSResidual(residualVec, t, dt, *_solution);
    } // for

    PYLITH_METHOD_END;
} // computeRHSResidual

// ----------------------------------------------------------------------
// Compute RHS Jacobian for G(t,s).
void
pylith::problems::Problem::computeRHSJacobian(PetscMat jacobianMat,
                                              PetscMat precondMat,
                                              const PylithReal t,
                                              const PylithReal dt,
                                              PetscVec solutionVec)
{ // computeRHSJacobian
    PYLITH_METHOD_BEGIN;
    PYLITH_JOURNAL_DEBUG("Problem::computeRHSJacobian(t="<<t<<", dt="<<dt<<", solutionVec="<<solutionVec<<", jacobianMat="<<jacobianMat<<", precondMat="<<precondMat<<")");

    // :KLUDGE: Should add check to see if we need to compute Jacobian

    // Update PyLith view of the solution.
    assert(_solution);
    _solution->scatterGlobalToLocal(solutionVec);

    // Sum Jacobian contributions across integrators.
    const size_t numIntegrators = _integrators.size();
    for (size_t i=0; i < numIntegrators; ++i) {
        _integrators[i]->computeRHSJacobian(jacobianMat, precondMat, t, dt, *_solution);
    } // for

    PYLITH_METHOD_END;
} // computeRHSJacobian


// ----------------------------------------------------------------------
// Compute LHS residual for F(t,s,\dot{s}).
void
pylith::problems::Problem::computeLHSResidual(PetscVec residualVec,
                                              const PylithReal t,
                                              const PylithReal dt,
                                              PetscVec solutionVec,
                                              PetscVec solutionDotVec)
{ // computeLHSResidual
    PYLITH_METHOD_BEGIN;
    PYLITH_JOURNAL_DEBUG("Problem::computeLHSResidual(t="<<t<<", dt="<<dt<<", solutionVec="<<solutionVec<<", solutionDotVec"<<solutionDotVec<<", residualVec="<<residualVec<<")");

    assert(residualVec);
    assert(solutionVec);
    assert(solutionDotVec);
    assert(_solution);

    // Update PyLith view of the solution.
    _solution->scatterGlobalToLocal(solutionVec);

    // Sum residual across integrators.
    const int numIntegrators = _integrators.size();
    assert(numIntegrators > 0); // must have at least 1 integrator
    for (int i=0; i < numIntegrators; ++i) {
        _integrators[i]->computeLHSResidual(residualVec, t, dt, *_solution, solutionDotVec);
    } // for

    PYLITH_METHOD_END;
} // computeLHSResidual


// ----------------------------------------------------------------------
// Compute LHS Jacobian for F(t,s,\dot{s}) for implicit time stepping.
void
pylith::problems::Problem::computeLHSJacobianImplicit(PetscMat jacobianMat,
                                                      PetscMat precondMat,
                                                      const PylithReal t,
                                                      const PylithReal dt,
                                                      const PylithReal tshift,
                                                      PetscVec solutionVec,
                                                      PetscVec solutionDotVec)
{ // computeLHSJacobianImplicit
    PYLITH_METHOD_BEGIN;
    PYLITH_JOURNAL_DEBUG("Problem::computeLHSJacobianImplicit(t="<<t<<", dt="<<dt<<", tshift="<<tshift<<", solutionVec="<<solutionVec<<", solutionDotVec="<<solutionDotVec<<", jacobianMat="<<jacobianMat<<", precondMat="<<precondMat<<")");

    // :KLUDGE: :TODO: Should add check to see if we need to compute Jacobian

    // Update PyLith view of the solution.
    assert(_solution);
    _solution->scatterGlobalToLocal(solutionVec);

    // Sum Jacobian contributions across integrators.
    const size_t numIntegrators = _integrators.size();
    for (size_t i=0; i < numIntegrators; ++i) {
        _integrators[i]->computeLHSJacobianImplicit(jacobianMat, precondMat, t, dt, tshift, *_solution, solutionDotVec);
    } // for

    PYLITH_METHOD_END;
} // computeLHSJacobianImplicit

// ----------------------------------------------------------------------
// Compute LHS Jacobian for F(t,s,\dot{s}) for explicit time stepping.
void
pylith::problems::Problem::computeLHSJacobianLumpedInv(const PylithReal t,
                                                       const PylithReal dt,
                                                       PetscVec solutionVec)
{ // computeLHSJacobianLumpedInv
    PYLITH_METHOD_BEGIN;
    PYLITH_JOURNAL_DEBUG("Problem::computeLHSJacobianLumpedInv(t="<<t<<", dt="<<dt<<", solutionVec="<<solutionVec<<")");

    // :KLUDGE: :TODO: Should add check to see if we need to compute Jacobian

    // Set jacobian to zero.
    assert(_jacobianLHSLumpedInv);
    _jacobianLHSLumpedInv->zeroAll();

    // Update PyLith view of the solution.
    assert(_solution);
    _solution->scatterGlobalToLocal(solutionVec);

    // Sum Jacobian contributions across integrators.
    const size_t numIntegrators = _integrators.size();
    for (size_t i=0; i < numIntegrators; ++i) {
        _integrators[i]->computeLHSJacobianLumpedInv(_jacobianLHSLumpedInv, t, dt, *_solution);
    } // for

    PYLITH_METHOD_END;
} // computeLHSJacobianLumpedInv

// End of file
