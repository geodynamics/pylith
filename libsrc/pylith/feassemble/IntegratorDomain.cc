// -*- C++ -*-
//
// ---------------------------------------------------------------------------------------------------------------------
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
// ---------------------------------------------------------------------------------------------------------------------
//

#include <portinfo>

#include "IntegratorDomain.hh" // implementation of object methods

#include "pylith/topology/Mesh.hh" // USES Mesh
#include "pylith/topology/Field.hh" // USES Field

#include "spatialdata/spatialdb/GravityField.hh" // HASA GravityField
#include "petscds.h" // USES PetscDS

#include "pylith/utils/journals.hh" // USES PYLITH_JOURNAL_*

#include <cassert> // USES assert()
#include <stdexcept> // USES std::runtime_error

extern "C" PetscErrorCode DMPlexComputeResidual_Internal(PetscDM dm,
                                                         IS cellIS,
                                                         PetscReal time,
                                                         PetscVec locX,
                                                         PetscVec locX_t,
                                                         PetscVec locF,
                                                         void *user);

extern "C" PetscErrorCode DMPlexComputeJacobian_Internal(PetscDM dm,
                                                         IS cellIS,
                                                         PetscReal t,
                                                         PetscReal X_tShift,
                                                         PetscVec X,
                                                         PetscVec X_t,
                                                         PetscMat Jac,
                                                         PetscMat JacP,
                                                         void *user);

// ---------------------------------------------------------------------------------------------------------------------
// Default constructor.
pylith::feassemble::IntegratorDomain::IntegratorDomain(const int dimension) :
    _dimension(dimension),
    _id(0) {}


// ---------------------------------------------------------------------------------------------------------------------
// Destructor.
pylith::feassemble::IntegratorDomain::~IntegratorDomain(void) {
    deallocate();
} // destructor


// ---------------------------------------------------------------------------------------------------------------------
// Deallocate PETSc and local data structures.
void
pylith::feassemble::IntegratorDomain::deallocate(void) {
    PYLITH_METHOD_BEGIN;

    pylith::feassemble::Integrator::deallocate();

    PYLITH_METHOD_END;
} // deallocate


// ---------------------------------------------------------------------------------------------------------------------
// Get spatial dimension of material.
int
pylith::feassemble::IntegratorDomain::getDimension(void) const {
    return _dimension;
} // getDimension


// ---------------------------------------------------------------------------------------------------------------------
// Set value of label material-id used to identify material cells.
void
pylith::feassemble::IntegratorDomain::setMaterialId(const int value) {
    PYLITH_JOURNAL_DEBUG("setMaterialId(value="<<value<<")");

    _materialId = value;
} // setMaterialId


// ---------------------------------------------------------------------------------------------------------------------
// Get value of label material-id used to identify material cells.
int
pylith::feassemble::IntegratorDomain::getMaterialId(void) const {
    return _id;
} // getMaterialId


// ---------------------------------------------------------------------------------------------------------------------
// Set gravity field.
void
pylith::feassemble::IntegratorDomain::setGravityField(spatialdata::spatialdb::GravityField* const g) {
    _gravityField = g;
} // setGravityField


// ---------------------------------------------------------------------------------------------------------------------
// Get mesh associated with integration domain.
const pylith::topology::Mesh&
pylith::feassemble::IntegratorDomain::getIntegrationDomainMesh(void) const {
    assert(_auxField);
    return _auxField->mesh();
} // domainMesh


// ---------------------------------------------------------------------------------------------------------------------
void
pylith::feassemble::IntegratorDomain::setKernelsRHSResidual(const std::vector<ResidualKernels>& kernels) {
    PYLITH_METHOD_BEGIN;
    PYLITH_JOURNAL_DEBUG("setKernelsRHSResidual(# kernels="<<kernels.size()<<")");

    _kernelsRHSResidual = kernels;

    PYLITH_METHOD_END;
} // setKernelsRHSResidual


// ---------------------------------------------------------------------------------------------------------------------
void
pylith::feassemble::IntegratorDomain::setKernelsRHSJacobianl(const std::vector<JacobianlKernels>& kernels) {
    PYLITH_METHOD_BEGIN;
    PYLITH_JOURNAL_DEBUG("setKernelsRHSJacobianl(# kernels="<<kernels.size()<<")");

    _kernelsRHSJacobianl = kernels;

    PYLITH_METHOD_END;
} // setKernelsRHSJacobianl


// ---------------------------------------------------------------------------------------------------------------------
void
pylith::feassemble::IntegratorDomain::setKernelsLHSResidual(const std::vector<ResidualKernels>& kernels) {
    PYLITH_METHOD_BEGIN;
    PYLITH_JOURNAL_DEBUG("setKernelsLHSResidual(# kernels="<<kernels.size()<<")");

    _kernelsLHSResidual = kernels;

    PYLITH_METHOD_END;
} // setKernelsLHSResidual


// ---------------------------------------------------------------------------------------------------------------------
void
pylith::feassemble::IntegratorDomain::setKernelsLHSJacobian(const std::vector<JacobianKernels>& kernels) {
    PYLITH_METHOD_BEGIN;
    PYLITH_JOURNAL_DEBUG("setKernelsLHSJacobian(# kernels="<<kernels.size()<<")");

    _kernelsLHSJacobian = kernels;

    PYLITH_METHOD_END;
} // setKernelsLHSJacobian


// ---------------------------------------------------------------------------------------------------------------------
void
pylith::feassemble::IntegratorDomain::setKernelsUpdateStateVars(const std::vector<ProjectKernels>& kernels) {
    PYLITH_METHOD_BEGIN;
    PYLITH_JOURNAL_DEBUG("setKernelsUpdateStateVars(# kernels="<<kernels.size()<<")");

    _kernelsUps = kernels;

    PYLITH_METHOD_END;
} // setKernelsUpdateStateVars


// ---------------------------------------------------------------------------------------------------------------------
void
pylith::feassemble::IntegratorDomain::setKernelsDerivedField(const std::vector<ProjectKernels>& kernels) {
    PYLITH_METHOD_BEGIN;
    PYLITH_JOURNAL_DEBUG("setKernelsDerivedField(# kernels="<<kernels.size()<<")");

    _kernelsDerivedField = kernels;

    PYLITH_METHOD_END;
} // setKernelsDerivedField


// ---------------------------------------------------------------------------------------------------------------------
// Compute RHS residual for G(t,s).
void
pylith::feassemble::IntegratorDomain::computeRHSResidual(pylith::topology::Field* residual,
                                                         const PylithReal t,
                                                         const PylithReal dt,
                                                         const pylith::topology::Field& solution) {
    PYLITH_METHOD_BEGIN;
    PYLITH_JOURNAL_DEBUG("computeRHSResidual(residual="<<residual<<", t="<<t<<", dt="<<dt<<", solution="<<solution.label()<<")");

    if (!_kernelsRHSResidual) { PYLITH_METHOD_END;}

    _setCurrentKernelsConstants(solution, dt);

    pylith::topology::Field solutionDot(solution.mesh()); // No dependence on time derivative of solution in RHS.
    solutionDot.label("solution_dot");
    _computeResidual(residual, _kernelsRHSResidual, t, dt, solution, solutionDot);

    PYLITH_METHOD_END;
} // computeRHSResidual


// ---------------------------------------------------------------------------------------------------------------------
// Compute RHS Jacobian for G(t,s).
void
pylith::feassemble::IntegratorDomain::computeRHSJacobian(PetscMat jacobianMat,
                                                         PetscMat precondMat,
                                                         const PylithReal t,
                                                         const PylithReal dt,
                                                         const pylith::topology::Field& solution) {
    PYLITH_METHOD_BEGIN;
    PYLITH_JOURNAL_DEBUG("computeRHSJacobian(jacobianMat="<<jacobianMat<<", precondMat="<<precondMat<<", t="<<t<<", dt="<<dt<<", solution="<<solution.label()<<")");

    if (!_kernelsRHSJacobian) { PYLITH_METHOD_END;}

    _setCurrentKernelsConstants(solution, dt);

    pylith::topology::Field solutionDot(solution.mesh()); // No dependence on time derivative of solution in RHS.
    solutionDot.label("solution_dot");
    const PylithReal s_tshift = 0.0; // No dependence on time derivative of solution in RHS, so shift isn't applicable.
    _computeJacobian(jacobianMat, precondMat, t_kernelsRHSJacobian, dt, s_tshift, solution, solutionDot);
    _needNewRHSJacobian = false;

    PYLITH_METHOD_END;
} // computeRHSJacobian


// ---------------------------------------------------------------------------------------------------------------------
// Compute LHS residual for F(t,s,\dot{s}).
void
pylith::feassemble::IntegratorDomain::computeLHSResidual(pylith::topology::Field* residual,
                                                         const PylithReal t,
                                                         const PylithReal dt,
                                                         const pylith::topology::Field& solution,
                                                         const pylith::topology::Field& solutionDot) {
    PYLITH_METHOD_BEGIN;
    PYLITH_JOURNAL_DEBUG("computeLHSResidual(residual="<<residual<<", t="<<t<<", dt="<<dt<<", solution="<<solution.label()<<")");

    if (!_kernelsLHSResidual) { PYLITH_METHOD_END;}

    _setCurrentKernelsConstants(solution, dt);

    _computeResidual(residual, _kernelsLHSResidualt, dt, solution, solutionDot);
} // computeLHSResidual


// ---------------------------------------------------------------------------------------------------------------------
// Compute LHS Jacobian for F(t,s,\dot{s}).
void
pylith::feassemble::IntegratorDomain::computeLHSJacobian(PetscMat jacobianMat,
                                                         PetscMat precondMat,
                                                         const PylithReal t,
                                                         const PylithReal dt,
                                                         const PylithReal s_tshift,
                                                         const pylith::topology::Field& solution,
                                                         const pylith::topology::Field& solutionDot) {
    PYLITH_METHOD_BEGIN;
    PYLITH_JOURNAL_DEBUG("computeLHSJacobian(jacobianMat="<<jacobianMat<<", precondMat="<<precondMat<<", t="<<t<<", dt="<<dt<<", solution="<<solution.label()<<", solutionDot="<<solutionDot.label()<<")");

    if (!_kernelsLHSJacobian) { PYLITH_METHOD_END;}

    _setCurrentKernelsJacobian(solution, _kernelsLHSJacobian);
    _setCurrentKernelsConstants(solution, dt);

    _computeJacobian(jacobianMat, precondMat, _kernelsLHSJacobian, t, dt, s_tshift, solution, solutionDot);
    _needNewLHSJacobian = false;

    PYLITH_METHOD_END;
} // computeLHSJacobian


// ---------------------------------------------------------------------------------------------------------------------
// Compute LHS Jacobian for F(t,s,\dot{s}).
void
pylith::feassemble::IntegratorDomain::computeLHSJacobianLumpedInv(pylith::topology::Field* jacobianInv,
                                                                  const PylithReal t,
                                                                  const PylithReal dt,
                                                                  const PylithReal s_tshift,
                                                                  const pylith::topology::Field& solution){
    PYLITH_METHOD_BEGIN;
    PYLITH_JOURNAL_DEBUG("computeLHSJacobianLumpedInv(jacobianInv="<<jacobianInv<<", t="<<t<<", dt="<<dt<<", solution="<<solution.label()<<")");

    assert(jacobianInv);

    _setCurrentKernelsConstants(solution, dt);

    PetscDS prob = NULL;
    PetscErrorCode err;

    PetscDM dmSoln = solution.dmMesh();
    PetscDM dmAux = _auxField->dmMesh();

    // Set pointwise function (kernels) in DS
    err = DMGetDS(dmSoln, &prob);PYLITH_CHECK_ERROR(err);
    const std::vector<JacobianKernels>& kernels = _kernelsLHSJacobian;
    for (size_t i = 0; i < kernels.size(); ++i) {
        const PetscInt i_fieldTrial = solution.subfieldInfo(kernels.fieldTrial).index;
        const PetscInt i_fieldBasis = solution.subfieldInfo(kernels.fieldBasis).index;
        err = PetscDSSetJacobian(prob, i_fieldTrial, i_fieldBasis, kernels.j0, kernels.j1, kernels.j2, kernels.j3);PYLITH_CHECK_ERROR(err);
    } // for

    // Get auxiliary data
    err = PetscObjectCompose((PetscObject) dmSoln, "dmAux", (PetscObject) dmAux);PYLITH_CHECK_ERROR(err);
    err = PetscObjectCompose((PetscObject) dmSoln, "A", (PetscObject) auxField()->localVector());PYLITH_CHECK_ERROR(err);

    PetscVec vecRowSum = NULL;
    err = DMGetGlobalVector(dmSoln, &vecRowSum);PYLITH_CHECK_ERROR(err);
    err = VecSet(vecRowSum, 1.0);PYLITH_CHECK_ERROR(err);

    // Compute the local Jacobian action
    PetscDMLabel dmLabel = NULL;
    PetscIS cells = NULL;
    PetscInt cStart = 0, cEnd = 0;
    err = DMGetLabel(dmSoln, "material-id", &dmLabel);PYLITH_CHECK_ERROR(err);
    err = DMLabelGetStratumBounds(dmLabel, id(), &cStart, &cEnd);PYLITH_CHECK_ERROR(err);
    err = ISCreateStride(PETSC_COMM_SELF, cEnd-cStart, cStart, 1, &cells);PYLITH_CHECK_ERROR(err);

    err = DMPlexComputeJacobianAction(dmSoln, cells, t, s_tshift, vecRowSum, NULL, vecRowSum, jacobianInv->localVector(), NULL);PYLITH_CHECK_ERROR(err);
    err = ISDestroy(&cells);PYLITH_CHECK_ERROR(err);

    // Compute the Jacobian inverse.
    err = VecReciprocal(jacobianInv->localVector());PYLITH_CHECK_ERROR(err);

    _needNewLHSJacobian = false;

    PYLITH_METHOD_END;
} // computeLHSJacobianLumpedInv


// ---------------------------------------------------------------------------------------------------------------------
// Compute residual using current kernels.
void
pylith::feassemble::IntegratorDomain::_computeResidual(pylith::topology::Field* residual,
                                                       const std::vector<ResidualKernels>& kernels,
                                                       const PylithReal t,
                                                       const PylithReal dt,
                                                       const pylith::topology::Field& solution,
                                                       const pylith::topology::Field& solutionDot) {
    PYLITH_METHOD_BEGIN;
    PYLITH_JOURNAL_DEBUG("_computeResidual(residual="<<residual<<", t="<<t<<", dt="<<dt<<", solution="<<solution.label()<<", solutionDot="<<solutionDot.label()<<")");

    assert(residual);
    assert(_auxField);

    PetscDS prob = NULL;
    PetscIS cells = NULL;
    PetscInt cStart = 0, cEnd = 0;
    PetscErrorCode err;

    PetscDM dmSoln = solution.dmMesh();
    PetscDM dmAux = _auxField->dmMesh();
    PetscDMLabel dmLabel = NULL;

    // Set pointwise function (kernels) in DS
    err = DMGetDS(dmSoln, &prob);PYLITH_CHECK_ERROR(err);
    for (size_t i = 0; i < kernels.size(); ++i) {
        const PetscInt i_field = solution.subfieldInfo(kernels.fieldTrial).index;
        err = PetscDSSetJacobian(prob, i_field, kernels.r0, kernels.r1);PYLITH_CHECK_ERROR(err);
    } // for

    // Get auxiliary data
    err = PetscObjectCompose((PetscObject) dmSoln, "dmAux", (PetscObject) dmAux);PYLITH_CHECK_ERROR(err);
    err = PetscObjectCompose((PetscObject) dmSoln, "A", (PetscObject) _auxField->localVector());PYLITH_CHECK_ERROR(err);

    // Compute the local residual
    assert(solution.localVector());
    assert(residual->localVector());
    err = DMGetLabel(dmSoln, "material-id", &dmLabel);PYLITH_CHECK_ERROR(err);
    err = DMLabelGetStratumBounds(dmLabel, id(), &cStart, &cEnd);PYLITH_CHECK_ERROR(err);
    assert(cEnd > cStart); // Double-check that this material has cells.

    PYLITH_JOURNAL_DEBUG("DMPlexComputeResidual_Internal() with material-id '"<<id()<<"' and cells ["<<cStart<<","<<cEnd<<".");
    err = ISCreateStride(PETSC_COMM_SELF, cEnd-cStart, cStart, 1, &cells);PYLITH_CHECK_ERROR(err);
    err = DMPlexComputeResidual_Internal(dmSoln, cells, PETSC_MIN_REAL, solution.localVector(), solutionDot.localVector(), residual->localVector(), NULL);PYLITH_CHECK_ERROR(err);
    err = ISDestroy(&cells);PYLITH_CHECK_ERROR(err);

    PYLITH_METHOD_END;
} // _computeResidual


// ---------------------------------------------------------------------------------------------------------------------
// Compute Jacobian using current kernels.
void
pylith::feassemble::IntegratorDomain::_computeJacobian(PetscMat jacobianMat,
                                                       PetscMat precondMat,
                                                       const std::vector<JacobianKernels>& kernels,
                                                       const PylithReal t,
                                                       const PylithReal dt,
                                                       const PylithReal s_tshift,
                                                       const pylith::topology::Field& solution,
                                                       const pylith::topology::Field& solutionDot) {
    PYLITH_METHOD_BEGIN;
    PYLITH_JOURNAL_DEBUG("_computeJacobian(jacobianMat="<<jacobianMat<<", precondMat="<<precondMat<<", t="<<t<<", dt="<<dt<<", s_tshift="<<s_tshift<<", solution="<<solution.label()<<", solutionDot="<<solutionDot.label()<<")");

    assert(jacobianMat);
    assert(precondMat);
    assert(_auxField);

    PetscDS prob = NULL;
    PetscIS cells = NULL;
    PetscInt cStart = 0, cEnd = 0;
    PetscErrorCode err;
    PetscDM dmMesh = solution.dmMesh();
    PetscDM dmAux = _auxField->dmMesh();
    PetscDMLabel dmLabel = NULL;

    // Set pointwise function (kernels) in DS
    err = DMGetDS(dmSoln, &prob);PYLITH_CHECK_ERROR(err);
    for (size_t i = 0; i < kernels.size(); ++i) {
        const PetscInt i_fieldTrial = solution.subfieldInfo(kernels.fieldTrial).index;
        const PetscInt i_fieldBasis = solution.subfieldInfo(kernels.fieldBasis).index;
        err = PetscDSSetJacobian(prob, i_fieldTrial, i_fieldBasis, kernels.j0, kernels.j1, kernels.j2, kernels.j3);PYLITH_CHECK_ERROR(err);
    } // for

    // Get auxiliary data
    err = PetscObjectCompose((PetscObject) dmMesh, "dmAux", (PetscObject) dmAux);PYLITH_CHECK_ERROR(err);
    err = PetscObjectCompose((PetscObject) dmMesh, "A", (PetscObject) auxField()->localVector());PYLITH_CHECK_ERROR(err);

    // Compute the local Jacobian
    assert(solution.localVector());
    err = DMGetLabel(dmMesh, "material-id", &dmLabel);PYLITH_CHECK_ERROR(err);
    err = DMLabelGetStratumBounds(dmLabel, id(), &cStart, &cEnd);PYLITH_CHECK_ERROR(err);

    PYLITH_JOURNAL_DEBUG("DMPlexComputeJacobian_Internal() with material-id '"<<id()<<"' and cells ["<<cStart<< ","<<cEnd<<".");
    err = ISCreateStride(PETSC_COMM_SELF, cEnd-cStart, cStart, 1, &cells);PYLITH_CHECK_ERROR(err);
    err = DMPlexComputeJacobian_Internal(dmMesh, cells, t, s_tshift, solution.localVector(), solutionDot.localVector(), jacobianMat, precondMat, NULL);PYLITH_CHECK_ERROR(err);
    err = ISDestroy(&cells);PYLITH_CHECK_ERROR(err);

    PYLITH_METHOD_END;
} // _computeJacobian


// End of file
