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

#include "pylith/feassemble/UpdateStateVars.hh" // HOLDSA UpdateStateVars
#include "pylith/topology/Mesh.hh" // USES Mesh
#include "pylith/topology/MeshOps.hh" // USES createSubdomainMesh()
#include "pylith/topology/Field.hh" // USES Field
#include "pylith/topology/CoordsVisitor.hh" // USES CoordsVisitor

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
pylith::feassemble::IntegratorDomain::IntegratorDomain(pylith::problems::Physics* const physics) :
    Integrator(physics),
    _materialId(0),
    _materialMesh(NULL),
    _updateState(NULL) {
    GenericComponent::setName("integratordomain");
} // constructor


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

    delete _materialMesh;_materialMesh = NULL;
    delete _updateState;_updateState = NULL;

    PYLITH_METHOD_END;
} // deallocate


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
    return _materialId;
} // getMaterialId


// ---------------------------------------------------------------------------------------------------------------------
// Get mesh associated with integration domain.
const pylith::topology::Mesh&
pylith::feassemble::IntegratorDomain::getPhysicsDomainMesh(void) const {
    assert(_materialMesh);
    return *_materialMesh;
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
pylith::feassemble::IntegratorDomain::setKernelsRHSJacobian(const std::vector<JacobianKernels>& kernels) {
    PYLITH_METHOD_BEGIN;
    PYLITH_JOURNAL_DEBUG("setKernelsRHSJacobian(# kernels="<<kernels.size()<<")");

    _kernelsRHSJacobian = kernels;

    PYLITH_METHOD_END;
} // setKernelsRHSJacobian


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

    _kernelsUpdateStateVars = kernels;

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
// Initialize integration domain, auxiliary field, and derived field. Update observers.
void
pylith::feassemble::IntegratorDomain::initialize(const pylith::topology::Field& solution) {
    PYLITH_METHOD_BEGIN;
    PYLITH_JOURNAL_DEBUG("intialize(solution="<<solution.getLabel()<<")");

    delete _materialMesh;
    _materialMesh = pylith::topology::MeshOps::createSubdomainMesh(solution.mesh(), "material-id", _materialId, ":UNKOWN:");
    pylith::topology::CoordsVisitor::optimizeClosure(_materialMesh->dmMesh());

    Integrator::initialize(solution);

    if (_kernelsUpdateStateVars.size() > 0) {
        delete _updateState;_updateState = new pylith::feassemble::UpdateStateVars;assert(_updateState);

        assert(_auxiliaryField);
        _updateState->initialize(*_auxiliaryField);
    } // if

    pythia::journal::debug_t debug(GenericComponent::getName());
    if (debug.state()) {
        PYLITH_JOURNAL_DEBUG("Viewing auxiliary field.");
        assert(_auxiliaryField);
        _auxiliaryField->view("Auxiliary field");
    } // if

    PYLITH_METHOD_END;
} // initialize


// ---------------------------------------------------------------------------------------------------------------------
// Compute RHS residual for G(t,s).
void
pylith::feassemble::IntegratorDomain::computeRHSResidual(pylith::topology::Field* residual,
                                                         const PylithReal t,
                                                         const PylithReal dt,
                                                         const pylith::topology::Field& solution) {
    PYLITH_METHOD_BEGIN;
    PYLITH_JOURNAL_DEBUG("computeRHSResidual(residual="<<residual<<", t="<<t<<", dt="<<dt<<", solution="<<solution.getLabel()<<")");

    if (0 == _kernelsRHSResidual.size()) { PYLITH_METHOD_END;}

    _setKernelConstants(solution, dt);

    pylith::topology::Field solutionDot(solution.mesh()); // No dependence on time derivative of solution in RHS.
    solutionDot.setLabel("solution_dot");
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
    PYLITH_JOURNAL_DEBUG("computeRHSJacobian(jacobianMat="<<jacobianMat<<", precondMat="<<precondMat<<", t="<<t<<", dt="<<dt<<", solution="<<solution.getLabel()<<")");

    _needNewRHSJacobian = false;
    if (0 == _kernelsRHSJacobian.size()) { PYLITH_METHOD_END;}

    _setKernelConstants(solution, dt);

    pylith::topology::Field solutionDot(solution.mesh()); // No dependence on time derivative of solution in RHS.
    solutionDot.setLabel("solution_dot");
    const PylithReal s_tshift = 0.0; // No dependence on time derivative of solution in RHS, so shift isn't applicable.
    _computeJacobian(jacobianMat, precondMat, _kernelsRHSJacobian, t, dt, s_tshift, solution, solutionDot);

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
    PYLITH_JOURNAL_DEBUG("computeLHSResidual(residual="<<residual<<", t="<<t<<", dt="<<dt<<", solution="<<solution.getLabel()<<")");

    if (0 == _kernelsLHSResidual.size()) { PYLITH_METHOD_END;}

    _setKernelConstants(solution, dt);
    _computeResidual(residual, _kernelsLHSResidual, t, dt, solution, solutionDot);

    PYLITH_METHOD_END;
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
    PYLITH_JOURNAL_DEBUG("computeLHSJacobian(jacobianMat="<<jacobianMat<<", precondMat="<<precondMat<<", t="<<t<<", dt="<<dt<<", solution="<<solution.getLabel()<<", solutionDot="<<solutionDot.getLabel()<<")");

    _needNewLHSJacobian = false;
    if (0 == _kernelsLHSJacobian.size()) { PYLITH_METHOD_END;}

    _setKernelConstants(solution, dt);
    _computeJacobian(jacobianMat, precondMat, _kernelsLHSJacobian, t, dt, s_tshift, solution, solutionDot);

    PYLITH_METHOD_END;
} // computeLHSJacobian


// ---------------------------------------------------------------------------------------------------------------------
// Compute LHS Jacobian for F(t,s,\dot{s}).
void
pylith::feassemble::IntegratorDomain::computeLHSJacobianLumpedInv(pylith::topology::Field* jacobianInv,
                                                                  const PylithReal t,
                                                                  const PylithReal dt,
                                                                  const PylithReal s_tshift,
                                                                  const pylith::topology::Field& solution) {
    PYLITH_METHOD_BEGIN;
    PYLITH_JOURNAL_DEBUG("computeLHSJacobianLumpedInv(jacobianInv="<<jacobianInv<<", t="<<t<<", dt="<<dt<<", solution="<<solution.getLabel()<<")");

    assert(jacobianInv);

    _setKernelConstants(solution, dt);

    PetscErrorCode err;

    // :KLUDGE: Potentially we may have multiple PetscDS objects. This assumes that the first one (with a NULL label) is
    // the correct one.
    PetscDS prob = NULL;
    PetscDM dmSoln = solution.dmMesh();assert(dmSoln);
    err = DMGetDS(dmSoln, &prob);PYLITH_CHECK_ERROR(err);assert(prob);

    // Set pointwise function (kernels) in DS
    const std::vector<JacobianKernels>& kernels = _kernelsLHSJacobian;
    for (size_t i = 0; i < kernels.size(); ++i) {
        const PetscInt i_fieldTrial = solution.subfieldInfo(kernels[i].subfieldTrial.c_str()).index;
        const PetscInt i_fieldBasis = solution.subfieldInfo(kernels[i].subfieldBasis.c_str()).index;
        err = PetscDSSetJacobian(prob, i_fieldTrial, i_fieldBasis, kernels[i].j0, kernels[i].j1, kernels[i].j2, kernels[i].j3);PYLITH_CHECK_ERROR(err);
    } // for

    // Get auxiliary data
    PetscDM dmAux = _auxiliaryField->dmMesh();assert(dmAux);
    err = PetscObjectCompose((PetscObject) dmSoln, "dmAux", (PetscObject) dmAux);PYLITH_CHECK_ERROR(err);
    err = PetscObjectCompose((PetscObject) dmSoln, "A", (PetscObject) _auxiliaryField->localVector());PYLITH_CHECK_ERROR(err);

    PetscVec vecRowSum = NULL;
    err = DMGetLocalVector(dmSoln, &vecRowSum);PYLITH_CHECK_ERROR(err);
    err = VecSet(vecRowSum, 1.0);PYLITH_CHECK_ERROR(err);

    // Compute the local Jacobian action
    PetscDMLabel dmLabel = NULL;
    PetscIS cells = NULL;
    PetscInt cStart = 0, cEnd = 0;
    err = DMGetLabel(dmSoln, "material-id", &dmLabel);PYLITH_CHECK_ERROR(err);
    err = DMLabelGetStratumBounds(dmLabel, _materialId, &cStart, &cEnd);PYLITH_CHECK_ERROR(err);
    err = ISCreateStride(PETSC_COMM_SELF, cEnd-cStart, cStart, 1, &cells);PYLITH_CHECK_ERROR(err);

    err = DMPlexComputeJacobianAction(dmSoln, cells, t, s_tshift, vecRowSum, NULL, vecRowSum, jacobianInv->localVector(), NULL);PYLITH_CHECK_ERROR(err);
    err = ISDestroy(&cells);PYLITH_CHECK_ERROR(err);

    // Compute the Jacobian inverse.
    err = VecReciprocal(jacobianInv->localVector());PYLITH_CHECK_ERROR(err);

    _needNewLHSJacobianLumped = false;

    PYLITH_METHOD_END;
} // computeLHSJacobianLumpedInv


// ---------------------------------------------------------------------------------------------------------------------
// Update state variables as needed.
void
pylith::feassemble::IntegratorDomain::_updateStateVars(const PylithReal t,
                                                       const PylithReal dt,
                                                       const pylith::topology::Field& solution) {
    PYLITH_METHOD_BEGIN;
    PYLITH_JOURNAL_DEBUG("_updateStateVars(t="<<t<<", dt="<<dt<<", solution="<<solution.getLabel()<<")");

    if (0 == _kernelsUpdateStateVars.size()) {
        PYLITH_METHOD_END;
    } // if

    assert(_updateState);
    assert(_auxiliaryField);
    _updateState->prepare(_auxiliaryField);
    _setKernelConstants(solution, dt);

    // We assume order of the update state variable kernels matches
    // the order of the correspoinding subfields in the auxiliary
    // field.
    const size_t numKernels = _kernelsUpdateStateVars.size();
    PetscPointFunc* kernelsStateVars = (numKernels > 0) ? new PetscPointFunc[numKernels] : NULL;
    for (size_t iKernel = 0; iKernel < numKernels; ++iKernel) {
        kernelsStateVars[iKernel] = _kernelsUpdateStateVars[iKernel].f;
    } // for

    PetscErrorCode err = 0;
    PetscDM stateVarsDM = _updateState->stateVarsDM();
    err = PetscObjectCompose((PetscObject) stateVarsDM, "dmAux", (PetscObject) _auxiliaryField->dmMesh());PYLITH_CHECK_ERROR(err);
    err = PetscObjectCompose((PetscObject) stateVarsDM, "A", (PetscObject) _auxiliaryField->localVector());PYLITH_CHECK_ERROR(err);
    err = DMProjectFieldLocal(stateVarsDM, t, solution.localVector(), kernelsStateVars, INSERT_VALUES,
                              _updateState->stateVarsLocalVector());PYLITH_CHECK_ERROR(err);
    _updateState->restore(_auxiliaryField);

    delete[] kernelsStateVars;kernelsStateVars = NULL;

    PYLITH_METHOD_END;
} // _updateStateVars


// ---------------------------------------------------------------------------------------------------------------------
// Compute field derived from solution and auxiliary field.
void
pylith::feassemble::IntegratorDomain::_computeDerivedField(const PylithReal t,
                                                           const PylithReal dt,
                                                           const pylith::topology::Field& solution) {
    PYLITH_METHOD_BEGIN;
    PYLITH_JOURNAL_DEBUG("_computeDerivedField(t="<<t<<", dt="<<dt<<", solution="<<solution.getLabel()<<")");

    if (!_derivedField) {
        PYLITH_METHOD_END;
    } // if

    assert(_derivedField);
    _setKernelConstants(solution, dt);

    const size_t numKernels = _kernelsDerivedField.size();
    PetscPointFunc* kernelsArray = (numKernels > 0) ? new PetscPointFunc[numKernels] : NULL;
    for (size_t iKernel = 0; iKernel < numKernels; ++iKernel) {
        const pylith::topology::Field::SubfieldInfo& sinfo = _derivedField->subfieldInfo(_kernelsDerivedField[iKernel].subfield.c_str());
        kernelsArray[sinfo.index] = _kernelsDerivedField[iKernel].f;
    } // for

    PetscErrorCode err = 0;

    PetscDM derivedDM = _derivedField->dmMesh();
    assert(_auxiliaryField);
    err = PetscObjectCompose((PetscObject) derivedDM, "dmAux", (PetscObject) _auxiliaryField->dmMesh());PYLITH_CHECK_ERROR(err);
    err = PetscObjectCompose((PetscObject) derivedDM, "A", (PetscObject) _auxiliaryField->localVector());PYLITH_CHECK_ERROR(err);

    err = DMProjectFieldLocal(derivedDM, t, solution.localVector(), kernelsArray, INSERT_VALUES, _derivedField->localVector());PYLITH_CHECK_ERROR(err);
    delete[] kernelsArray;kernelsArray = NULL;

    pythia::journal::debug_t debug(GenericComponent::getName());
    if (debug.state()) {
        PYLITH_JOURNAL_DEBUG("Viewing derived field.");
        _derivedField->view("Derived field");
    } // if

    PYLITH_METHOD_END;
} // _computeDerivedField


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
    PYLITH_JOURNAL_DEBUG("_computeResidual(residual="<<residual<<", # kernels="<<kernels.size()<<", t="<<t<<", dt="<<dt<<", solution="<<solution.getLabel()<<", solutionDot="<<solutionDot.getLabel()<<")");

    assert(residual);
    assert(_auxiliaryField);

    PetscDS prob = NULL;
    PetscIS cells = NULL;
    PetscInt cStart = 0, cEnd = 0;
    PetscErrorCode err;

    PetscDM dmSoln = solution.dmMesh();
    PetscDM dmAux = _auxiliaryField->dmMesh();
    PetscDMLabel dmLabel = NULL;

    // :KLUDGE: Potentially we may have multiple PetscDS objects. This assumes that the first one (with a NULL label) is
    // the correct one.

    // Set pointwise function (kernels) in DS
    err = DMGetDS(dmSoln, &prob);PYLITH_CHECK_ERROR(err);
    for (size_t i = 0; i < kernels.size(); ++i) {
        const PetscInt i_field = solution.subfieldInfo(kernels[i].subfield.c_str()).index;
        err = PetscDSSetResidual(prob, i_field, kernels[i].r0, kernels[i].r1);PYLITH_CHECK_ERROR(err);
    } // for

    // Get auxiliary data
    err = PetscObjectCompose((PetscObject) dmSoln, "dmAux", (PetscObject) dmAux);PYLITH_CHECK_ERROR(err);
    err = PetscObjectCompose((PetscObject) dmSoln, "A", (PetscObject) _auxiliaryField->localVector());PYLITH_CHECK_ERROR(err);

    // Compute the local residual
    assert(solution.localVector());
    assert(residual->localVector());
    err = DMGetLabel(dmSoln, "material-id", &dmLabel);PYLITH_CHECK_ERROR(err);
    err = DMLabelGetStratumBounds(dmLabel, _materialId, &cStart, &cEnd);PYLITH_CHECK_ERROR(err);
    assert(cEnd > cStart); // Double-check that this material has cells.

    PYLITH_JOURNAL_DEBUG("DMPlexComputeResidual_Internal() with material-id '"<<_materialId<<"' and cells ["<<cStart<<","<<cEnd<<").");
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
    PYLITH_JOURNAL_DEBUG("_computeJacobian(jacobianMat="<<jacobianMat<<", precondMat="<<precondMat<<", # kernels="<<kernels.size()<<", t="<<t<<", dt="<<dt<<", s_tshift="<<s_tshift<<", solution="<<solution.getLabel()<<", solutionDot="<<solutionDot.getLabel()<<")");

    assert(jacobianMat);
    assert(precondMat);
    assert(_auxiliaryField);

    PetscDS prob = NULL;
    PetscIS cells = NULL;
    PetscInt cStart = 0, cEnd = 0;
    PetscErrorCode err;
    PetscDM dmMesh = solution.dmMesh();
    PetscDM dmAux = _auxiliaryField->dmMesh();
    PetscDMLabel dmLabel = NULL;

    // :KLUDGE: Potentially we may have multiple PetscDS objects. This assumes that the first one (with a NULL label) is
    // the correct one.

    // Set pointwise function (kernels) in DS
    err = DMGetDS(solution.dmMesh(), &prob);PYLITH_CHECK_ERROR(err);
    for (size_t i = 0; i < kernels.size(); ++i) {
        const PetscInt i_fieldTrial = solution.subfieldInfo(kernels[i].subfieldTrial.c_str()).index;
        const PetscInt i_fieldBasis = solution.subfieldInfo(kernels[i].subfieldBasis.c_str()).index;
        err = PetscDSSetJacobian(prob, i_fieldTrial, i_fieldBasis, kernels[i].j0, kernels[i].j1, kernels[i].j2, kernels[i].j3);PYLITH_CHECK_ERROR(err);
    } // for

    // Get auxiliary data
    err = PetscObjectCompose((PetscObject) dmMesh, "dmAux", (PetscObject) dmAux);PYLITH_CHECK_ERROR(err);
    err = PetscObjectCompose((PetscObject) dmMesh, "A", (PetscObject) _auxiliaryField->localVector());PYLITH_CHECK_ERROR(err);

    // Compute the local Jacobian
    assert(solution.localVector());
    err = DMGetLabel(dmMesh, "material-id", &dmLabel);PYLITH_CHECK_ERROR(err);
    err = DMLabelGetStratumBounds(dmLabel, _materialId, &cStart, &cEnd);PYLITH_CHECK_ERROR(err);

    PYLITH_JOURNAL_DEBUG("DMPlexComputeJacobian_Internal() with material-id '"<<_materialId<<"' and cells ["<<cStart<< ","<<cEnd<<".");
    err = ISCreateStride(PETSC_COMM_SELF, cEnd-cStart, cStart, 1, &cells);PYLITH_CHECK_ERROR(err);
    err = DMPlexComputeJacobian_Internal(dmMesh, cells, t, s_tshift, solution.localVector(), solutionDot.localVector(), jacobianMat, precondMat, NULL);PYLITH_CHECK_ERROR(err);
    err = ISDestroy(&cells);PYLITH_CHECK_ERROR(err);

    PYLITH_METHOD_END;
} // _computeJacobian


// End of file
