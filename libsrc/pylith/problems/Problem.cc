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

#include "pylith/materials/Material.hh" // USES Material
#include "pylith/faults/FaultCohesive.hh" // USES FaultCohesive
#include "pylith/bc/BoundaryCondition.hh" // USES BoundaryCondition
#include "pylith/feassemble/Integrator.hh" // USES Integrator
#include "pylith/feassemble/Constraint.hh" // USES Constraint
#include "pylith/feassemble/Observers.hh" // USES Observers
#include "pylith/topology/MeshOps.hh" // USES MeshOps
#include "pylith/topology/CoordsVisitor.hh" // USES CoordsVisitor::optimizeClosure()

#include "spatialdata/units/Nondimensional.hh" // USES Nondimensional
#include "spatialdata/spatialdb/GravityField.hh" // USES GravityField

#include "pylith/utils/error.hh" // USES PYLITH_CHECK_ERROR
#include "pylith/utils/journals.hh" // USES PYLITH_COMPONENT_*

#include <cassert> // USES assert()
#include <typeinfo> \
    // USES typeid()

// ----------------------------------------------------------------------
// Constructor
pylith::problems::Problem::Problem() :
    _solution(NULL),
    _solutionDot(NULL),
    _residual(NULL),
    _jacobianLHSLumpedInv(NULL),
    _normalizer(NULL),
    _gravityField(NULL),
    _observers(new pylith::feassemble::Observers),
    _solverType(LINEAR) { // constructor
} // constructor


// ---------------------------------------------------------------------------------------------------------------------
// Destructor
pylith::problems::Problem::~Problem(void) {
    deallocate();
} // destructor


// ---------------------------------------------------------------------------------------------------------------------
// Deallocate PETSc and local data structures.
void
pylith::problems::Problem::deallocate(void) {
    PYLITH_METHOD_BEGIN;

    _solution = NULL; // Held by Python. :KLUDGE: :TODO: Use shared pointer.
    delete _solutionDot;_solutionDot = NULL;
    delete _residual;_residual = NULL;
    delete _jacobianLHSLumpedInv;_jacobianLHSLumpedInv = NULL;
    delete _normalizer;_normalizer = NULL;
    _gravityField = NULL; // Held by Python. :KLUDGE: :TODO: Use shared pointer.

    PYLITH_METHOD_END;
} // deallocate


// ---------------------------------------------------------------------------------------------------------------------
// Set problem type.
void
pylith::problems::Problem::setSolverType(const SolverTypeEnum value) {
    PYLITH_COMPONENT_DEBUG("Problem::setSolverType(value="<<value<<")");

    _solverType = value;
} // setSolverType


// ---------------------------------------------------------------------------------------------------------------------
// Get problem type.
pylith::problems::Problem::SolverTypeEnum
pylith::problems::Problem::getSolverType(void) const {
    return _solverType;
} // getSolverType


// ---------------------------------------------------------------------------------------------------------------------
// Set manager of scales used to nondimensionalize problem.
void
pylith::problems::Problem::setNormalizer(const spatialdata::units::Nondimensional& dim) {
    PYLITH_COMPONENT_DEBUG("Problem::setNormalizer(dim="<<typeid(dim).name()<<")");

    if (!_normalizer) {
        _normalizer = new spatialdata::units::Nondimensional(dim);
    } else {
        *_normalizer = dim;
    } // if/else
} // setNormalizer


// ---------------------------------------------------------------------------------------------------------------------
// Set gravity field.
void
pylith::problems::Problem::setGravityField(spatialdata::spatialdb::GravityField* const g) {
    PYLITH_COMPONENT_DEBUG("Problem::setGravityField(g="<<typeid(*g).name()<<")");

    _gravityField = g;
} // setGravityField


// ----------------------------------------------------------------------
// Register observer to receive notifications.
void
pylith::problems::Problem::registerObserver(pylith::feassemble::Observer* observer) {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("registerObserver(observer="<<typeid(observer).name()<<")");

    assert(_observers);
    _observers->registerObserver(observer);

    PYLITH_METHOD_END;
} // registerObserver


// ----------------------------------------------------------------------
// Remove observer from receiving notifications.
void
pylith::problems::Problem::removeObserver(pylith::feassemble::Observer* observer) {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("removeObserver(observer="<<typeid(observer).name()<<")");

    assert(_observers);
    _observers->removeObserver(observer);

    PYLITH_METHOD_END;
} // removeObserver


// ---------------------------------------------------------------------------------------------------------------------
// Set solution field.
void
pylith::problems::Problem::setSolution(pylith::topology::Field* field) {
    PYLITH_COMPONENT_DEBUG("Problem::setSolution(field="<<typeid(*field).name()<<")");

    _solution = field;
} // setSolution


// ---------------------------------------------------------------------------------------------------------------------
// Set materials.
void
pylith::problems::Problem::setMaterials(pylith::materials::Material* materials[],
                                        const int numMaterials) {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("Problem::setMaterials("<<materials<<", numMaterials="<<numMaterials<<")");

    assert( (!materials && 0 == numMaterials) || (materials && 0 < numMaterials) );

    _materials.resize(numMaterials);
    for (int i = 0; i < numMaterials; ++i) {
        _materials[i] = materials[i];
    } // for

    PYLITH_METHOD_END;
} // setMaterials


// ---------------------------------------------------------------------------------------------------------------------
// Set materials.
void
pylith::problems::Problem::setBoundaryConditions(pylith::bc::BoundaryCondition* bc[],
                                                 const int numBC) {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("Problem::setBoundaryConditions("<<bc<<", numBC="<<numBC<<")");

    assert( (!bc && 0 == numBC) || (bc && 0 < numBC) );

    _bc.resize(numBC);
    for (int i = 0; i < numBC; ++i) {
        _bc[i] = bc[i];
    } // for

    PYLITH_METHOD_END;
} // setBoundaryConditions


// ---------------------------------------------------------------------------------------------------------------------
// Set materials.
void
pylith::problems::Problem::setInterfaces(pylith::faults::FaultCohesive* interfaces[],
                                         const int numInterfaces) {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("Problem::setInterfaces("<<interfaces<<", numInterfaces="<<numInterfaces<<")");

    assert( (!interfaces && 0 == numInterfaces) || (interfaces && 0 < numInterfaces) );

    _interfaces.resize(numInterfaces);
    for (int i = 0; i < numInterfaces; ++i) {
        _interfaces[i] = interfaces[i];
    } // for

    PYLITH_METHOD_END;
} // setInterfaces


// ----------------------------------------------------------------------
// Do minimal initialization.
void
pylith::problems::Problem::preinitialize(const pylith::topology::Mesh& mesh) {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("Problem::preinitialzie(mesh="<<typeid(mesh).name()<<")");

    assert(_normalizer);

    const size_t numMaterials = _materials.size();
    for (size_t i = 0; i < numMaterials; ++i) {
        assert(_materials[i]);
        _materials[i]->setNormalizer(*_normalizer);
        _materials[i]->setGravityField(_gravityField);
    } // for

    const size_t numInterfaces = _interfaces.size();
    for (size_t i = 0; i < numInterfaces; ++i) {
        assert(_interfaces[i]);
        _interfaces[i]->setNormalizer(*_normalizer);
    } // for

    const size_t numBC = _bc.size();
    for (size_t i = 0; i < numBC; ++i) {
        assert(_bc[i]);
        _bc[i]->setNormalizer(*_normalizer);
    } // for

    PYLITH_METHOD_END;
} // preinitialize


// ----------------------------------------------------------------------
// Verify configuration.
void
pylith::problems::Problem::verifyConfiguration(void) const {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("Problem::verifyConfiguration(void)");

    assert(_solution);

    // Check to make sure materials are compatible with the solution.
    const size_t numMaterials = _materials.size();
    for (size_t i = 0; i < numMaterials; ++i) {
        assert(_materials[i]);
        _materials[i]->verifyConfiguration(*_solution);
    } // for

    // Check to make sure interfaces are compatible with the solution.
    const size_t numInterfaces = _interfaces.size();
    for (size_t i = 0; i < numInterfaces; ++i) {
        assert(_interfaces[i]);
        _interfaces[i]->verifyConfiguration(*_solution);
    } // for

    // Check to make sure boundary conditions are compatible with the solution.
    const size_t numBC = _bc.size();
    for (size_t i = 0; i < numBC; ++i) {
        assert(_bc[i]);
        _bc[i]->verifyConfiguration(*_solution);
    } // for

    assert(_observers);
    _observers->verifyObservers(*_solution);

    PYLITH_METHOD_END;
} // verifyConfiguration


// ----------------------------------------------------------------------
// Initialize.
void
pylith::problems::Problem::initialize(void) {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("Problem::initialize()");

    assert(_solution);

    _createIntegrators();
    _createConstraints();

    const pylith::topology::Mesh& mesh = _solution->mesh();
    PetscDM dmMesh = mesh.dmMesh();assert(dmMesh);
    pylith::topology::CoordsVisitor::optimizeClosure(dmMesh);

    // Initialize integrators.
    const size_t numIntegrators = _integrators.size();
    for (size_t i = 0; i < numIntegrators; ++i) {
        assert(_integrators[i]);
        _integrators[i]->initialize(*_solution);
    } // for

    // Initialize constraints.
    const size_t numConstraints = _constraints.size();
    for (size_t i = 0; i < numConstraints; ++i) {
        assert(_constraints[i]);
        _constraints[i]->initialize(*_solution);
    } // for

    // Initialize solution field.
    if (_solution->hasSubfield("lagrange_multiplier_fault")) {
        _setupLagrangeMultiplier(_solution);
    } // if
    _solution->allocate();
    _solution->zeroLocal();
    _solution->createScatter(_solution->mesh(), "global");

    // Initialize residual.
    delete _residual;_residual = new pylith::topology::Field(_solution->mesh());assert(_residual);
    _residual->cloneSection(*_solution);
    _residual->label("residual");
    _solution->zeroLocal();

    PYLITH_METHOD_END;
} // initialize


// ----------------------------------------------------------------------
// Set solution values according to constraints (Dirichlet BC).
void
pylith::problems::Problem::setSolutionLocal(const PylithReal t,
                                            PetscVec solutionVec,
                                            PetscVec solutionDotVec) {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("setSolutionLocal(t="<<t<<", solutionVec="<<solutionVec<<", solutionDotVec="<<solutionDotVec<<")");

    // Update PyLith view of the solution.
    assert(_solution);
    _solution->scatterVectorToLocal(solutionVec);

    if (solutionDotVec) {
        if (!_solutionDot) {
            _solutionDot = new pylith::topology::Field(_solution->mesh());
            _solutionDot->cloneSection(*_solution);
            _solutionDot->label("solutionDot");
        } // if
        _solutionDot->scatterVectorToLocal(solutionDotVec);
    } // if

    const size_t numConstraints = _constraints.size();
    for (size_t i = 0; i < numConstraints; ++i) {
        _constraints[i]->setSolution(_solution, t);
    } // for

    // _solution->view("SOLUTION AFTER SETTING VALUES");

    PYLITH_METHOD_END;
} // setSolutionLocal


// ----------------------------------------------------------------------
// Compute RHS residual for G(t,s).
void
pylith::problems::Problem::computeRHSResidual(PetscVec residualVec,
                                              const PylithReal t,
                                              const PylithReal dt,
                                              PetscVec solutionVec) {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("Problem::computeRHSResidual(t="<<t<<", dt="<<dt<<", solutionVec="<<solutionVec<<", residualVec="<<residualVec<<")");

    assert(residualVec);
    assert(solutionVec);
    assert(_solution);

    // Update PyLith view of the solution.
    PetscVec solutionDotVec = NULL;
    setSolutionLocal(t, solutionVec, solutionDotVec);

    // Sum residual contributions across integrators.
    _residual->zeroLocal();
    const size_t numIntegrators = _integrators.size();
    assert(numIntegrators > 0); // must have at least 1 integrator
    for (size_t i = 0; i < numIntegrators; ++i) {
        _integrators[i]->computeRHSResidual(_residual, t, dt, *_solution);
    } // for

    // Assemble residual values across processes.
    PetscErrorCode err = VecSet(residualVec, 0.0);PYLITH_CHECK_ERROR(err); // Move to TSComputeIFunction()?
    _residual->scatterLocalToVector(residualVec, ADD_VALUES);

    PYLITH_METHOD_END;
} // computeRHSResidual


// ----------------------------------------------------------------------
// Compute RHS Jacobian for G(t,s).
void
pylith::problems::Problem::computeRHSJacobian(PetscMat jacobianMat,
                                              PetscMat precondMat,
                                              const PylithReal t,
                                              const PylithReal dt,
                                              PetscVec solutionVec) {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("Problem::computeRHSJacobian(t="<<t<<", dt="<<dt<<", solutionVec="<<solutionVec<<", jacobianMat="<<jacobianMat<<", precondMat="<<precondMat<<")");

    assert(jacobianMat);
    assert(precondMat);
    assert(solutionVec);
    assert(_solution);

    const size_t numIntegrators = _integrators.size();

    PetscErrorCode err;
    err = MatZeroEntries(precondMat);PYLITH_CHECK_ERROR(err);

    // Check to see if we need to compute RHS Jacobian.
    bool needNewRHSJacobian = false;
    for (size_t i = 0; i < numIntegrators; ++i) {
        if (_integrators[i]->needNewRHSJacobian()) {
            needNewRHSJacobian = true;
            break;
        } // if
    } // for
    if (!needNewRHSJacobian) {
        PYLITH_METHOD_END;
    } // if

    // Update PyLith view of the solution.
    const PetscVec solutionDotVec = NULL;
    setSolutionLocal(t, solutionVec, solutionDotVec);

    // Sum Jacobian contributions across integrators.
    for (size_t i = 0; i < numIntegrators; ++i) {
        _integrators[i]->computeRHSJacobian(jacobianMat, precondMat, t, dt, *_solution);
    } // for

    // Solver handles assembly.

    PYLITH_METHOD_END;
} // computeRHSJacobian


// ----------------------------------------------------------------------
// Compute LHS residual for F(t,s,\dot{s}).
void
pylith::problems::Problem::computeLHSResidual(PetscVec residualVec,
                                              const PylithReal t,
                                              const PylithReal dt,
                                              PetscVec solutionVec,
                                              PetscVec solutionDotVec) {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("Problem::computeLHSResidual(t="<<t<<", dt="<<dt<<", solutionVec="<<solutionVec<<", solutionDotVec="<<solutionDotVec<<", residualVec="<<residualVec<<")");

    assert(residualVec);
    assert(solutionVec);
    assert(solutionDotVec);
    assert(_solution);

    // Update PyLith view of the solution.
    setSolutionLocal(t, solutionVec, solutionDotVec);

    // Sum residual across integrators.
    _residual->zeroLocal();
    const int numIntegrators = _integrators.size();
    assert(numIntegrators > 0); // must have at least 1 integrator
    for (int i = 0; i < numIntegrators; ++i) {
        _integrators[i]->computeLHSResidual(_residual, t, dt, *_solution, *_solutionDot);
    } // for

    // Assemble residual values across processes.
    PetscErrorCode err = VecSet(residualVec, 0.0);PYLITH_CHECK_ERROR(err); // Move to TSComputeIFunction()?
    _residual->scatterLocalToVector(residualVec, ADD_VALUES);

    PYLITH_METHOD_END;
} // computeLHSResidual


// ----------------------------------------------------------------------
// Compute LHS Jacobian for F(t,s,\dot{s}) for implicit time stepping.
void
pylith::problems::Problem::computeLHSJacobian(PetscMat jacobianMat,
                                              PetscMat precondMat,
                                              const PylithReal t,
                                              const PylithReal dt,
                                              const PylithReal s_tshift,
                                              PetscVec solutionVec,
                                              PetscVec solutionDotVec) {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("Problem::computeLHSJacobian(t="<<t<<", dt="<<dt<<", s_tshift="<<s_tshift<<", solutionVec="<<solutionVec<<", solutionDotVec="<<solutionDotVec<<", jacobianMat="<<jacobianMat<<", precondMat="<<precondMat<<")");

    assert(jacobianMat);
    assert(precondMat);
    assert(solutionVec);
    assert(solutionDotVec);
    assert(s_tshift > 0);

    const size_t numIntegrators = _integrators.size();

    // Check to see if we need to compute RHS Jacobian.
    bool needNewLHSJacobian = false;
    for (size_t i = 0; i < numIntegrators; ++i) {
        if (_integrators[i]->needNewLHSJacobian()) {
            needNewLHSJacobian = true;
            break;
        } // if
    } // for
    if (!needNewLHSJacobian) {
        PYLITH_METHOD_END;
    } // if

    // Update PyLith view of the solution.
    setSolutionLocal(t, solutionVec, solutionDotVec);

    // Sum Jacobian contributions across integrators.
    for (size_t i = 0; i < numIntegrators; ++i) {
        _integrators[i]->computeLHSJacobian(jacobianMat, precondMat, t, dt, s_tshift, *_solution, *_solutionDot);
    } // for

    // Solver handles assembly.

    PYLITH_METHOD_END;
} // computeLHSJacobianImplicit


// ----------------------------------------------------------------------
// Compute inverse of LHS Jacobian for F(t,s,\dot{s}) for explicit time stepping.
void
pylith::problems::Problem::computeLHSJacobianLumpedInv(const PylithReal t,
                                                       const PylithReal dt,
                                                       const PylithReal s_tshift,
                                                       PetscVec solutionVec) {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("Problem::computeLHSJacobianLumpedInv(t="<<t<<", dt="<<dt<<", s_tshift="<<s_tshift<<", solutionVec="<<solutionVec<<")");

    assert(solutionVec);
    assert(_solution);
    assert(_jacobianLHSLumpedInv);
    assert(s_tshift > 0);

    const size_t numIntegrators = _integrators.size();

    // Check to see if we need to compute LHS Jacobian.
    bool needNewLHSJacobian = false;
    for (size_t i = 0; i < numIntegrators; ++i) {
        if (_integrators[i]->needNewLHSJacobian()) {
            needNewLHSJacobian = true;
            break;
        } // if
    } // for
    if (!needNewLHSJacobian) {
        PYLITH_METHOD_END;
    } // if

    // Set jacobian to zero.
    _jacobianLHSLumpedInv->zeroLocal();

    // Update PyLith view of the solution.
    const PetscVec solutionDotVec = NULL;
    setSolutionLocal(t, solutionVec, solutionDotVec);

    // Sum Jacobian contributions across integrators.
    for (size_t i = 0; i < numIntegrators; ++i) {
        _integrators[i]->computeLHSJacobianLumpedInv(_jacobianLHSLumpedInv, t, dt, s_tshift, *_solution);
    } // for

    // No need to assemble inverse of lumped Jacobian across processes, because it contributes to residual.

    PYLITH_METHOD_END;
} // computeLHSJacobianLumpedInv


// ---------------------------------------------------------------------------------------------------------------------
// Check material and interface ids.
void
pylith::problems::Problem::_checkMaterialIds(void) {
    const size_t numMaterials = _materials.size();
    const size_t numInterfaces = _interfaces.size();

    pylith::int_array materialIds(numMaterials + numInterfaces);
    size_t count = 0;
    for (size_t i = 0; i < numMaterials; ++i) {
        assert(_materials[i]);
        materialIds[count++] = _materials[i]->getMaterialId();
    } // for
    for (size_t i = 0; i < numInterfaces; ++i) {
        assert(_interfaces[i]);
        materialIds[count++] = _interfaces[i]->getInterfaceId();
    } // for

    pylith::topology::MeshOps::checkMaterialIds(_solution->mesh(), materialIds);
} // _checkMaterialIds


// ---------------------------------------------------------------------------------------------------------------------
// Create array of integrators from materials, interfaces, and boundary conditions.
void
pylith::problems::Problem::_createIntegrators(void) {
    PYLITH_METHOD_BEGIN;

    const size_t numMaterials = _materials.size();
    const size_t numInterfaces = _interfaces.size();
    const size_t numBC = _bc.size();

    const size_t maxSize = numMaterials + numInterfaces + numBC;
    _integrators.resize(maxSize);
    size_t count = 0;

    for (size_t i = 0; i < numMaterials; ++i) {
        assert(_materials[i]);
        pylith::feassemble::Integrator* integrator = _materials[i]->createIntegrator(*_solution);
        assert(count < maxSize);
        if (integrator) { _integrators[count++] = integrator;}
    } // for

    for (size_t i = 0; i < numInterfaces; ++i) {
        assert(_interfaces[i]);
        pylith::feassemble::Integrator* integrator = _interfaces[i]->createIntegrator(*_solution);
        assert(count < maxSize);
        if (integrator) { _integrators[count++] = integrator;}
    } // for

    // Check to make sure boundary conditions are compatible with the solution.
    for (size_t i = 0; i < numBC; ++i) {
        assert(_bc[i]);
        pylith::feassemble::Integrator* integrator = _bc[i]->createIntegrator(*_solution);
        assert(count < maxSize);
        if (integrator) { _integrators[count++] = integrator;}
    } // for

    _integrators.resize(count);

    PYLITH_METHOD_END;
} // _createIntegrators


// ---------------------------------------------------------------------------------------------------------------------
// Create array of constraints from materials, interfaces, and boundary conditions.
void
pylith::problems::Problem::_createConstraints(void) {
    PYLITH_METHOD_BEGIN;

    const size_t numMaterials = _materials.size();
    const size_t numInterfaces = _interfaces.size();
    const size_t numBC = _bc.size();

    const size_t maxSize = numMaterials + numInterfaces + numBC;
    _constraints.resize(maxSize);
    size_t count = 0;

    for (size_t i = 0; i < numMaterials; ++i) {
        assert(_materials[i]);
        pylith::feassemble::Constraint* constraint = _materials[i]->createConstraint(*_solution);
        assert(count < maxSize);
        if (constraint) { _constraints[count++] = constraint;}
    } // for

    for (size_t i = 0; i < numInterfaces; ++i) {
        assert(_interfaces[i]);
        pylith::feassemble::Constraint* constraint = _interfaces[i]->createConstraint(*_solution);
        assert(count < maxSize);
        if (constraint) { _constraints[count++] = constraint;}
    } // for

    for (size_t i = 0; i < numBC; ++i) {
        assert(_bc[i]);
        pylith::feassemble::Constraint* constraint = _bc[i]->createConstraint(*_solution);
        assert(count < maxSize);
        if (constraint) { _constraints[count++] = constraint;}
    } // for

    _constraints.resize(count);

    PYLITH_METHOD_END;
} // _createConstraints


// ---------------------------------------------------------------------------------------------------------------------
// Setup field so Lagrange multiplier subfield is limited to degrees of freedom associated with the cohesive cells.
void
pylith::problems::Problem::_setupLagrangeMultiplier(pylith::topology::Field* solution) {
    PYLITH_METHOD_BEGIN;

    assert(_solution->hasSubfield("lagrange_multiplier_fault"));

    PetscDMLabel cohesiveLabel = NULL;
    PylithInt dim = 0;
    PylithInt pStart = 0;
    PylithInt pEnd = 0;
    PetscErrorCode err;

    PetscDM dmSoln = _solution->dmMesh();assert(dmSoln);
    err = DMGetDimension(dmSoln, &dim);PYLITH_CHECK_ERROR(err);
    PylithInt* pMax = dim+1 > 0 ? new PylithInt[dim+1] : NULL;
    err = DMPlexGetHybridBounds(dmSoln, dim > 0 ? &pMax[dim] : NULL, dim > 1 ? &pMax[dim-1] : NULL, dim > 2 ? &pMax[1] : NULL, &pMax[0]);PYLITH_CHECK_ERROR(err);
    err = DMCreateLabel(dmSoln, "cohesive interface");PYLITH_CHECK_ERROR(err);
    err = DMGetLabel(dmSoln, "cohesive interface", &cohesiveLabel);PYLITH_CHECK_ERROR(err);
    for (PylithInt iDim = 0; iDim <= dim; ++iDim) {
        err = DMPlexGetDepthStratum(dmSoln, iDim, &pStart, &pEnd);PYLITH_CHECK_ERROR(err);
        pStart = pMax[iDim] < 0 ? pEnd : pMax[iDim];
        for (PylithInt p = pStart; p < pEnd; ++p) {
            err = DMLabelSetValue(cohesiveLabel, p, 1);PYLITH_CHECK_ERROR(err);
        } // for
    } // for
    delete[] pMax;pMax = NULL;

    PetscDS prob = NULL;
    err = DMGetDS(dmSoln, &prob);
    err = PetscDSSetLabel(prob, 1, cohesiveLabel);PYLITH_CHECK_ERROR(err);

    PYLITH_METHOD_END;
} // _setupLagrangeMultiplier


// End of file
