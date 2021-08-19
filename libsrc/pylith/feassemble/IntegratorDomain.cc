// -*- C++ -*-
//
// ---------------------------------------------------------------------------------------------------------------------
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
// ---------------------------------------------------------------------------------------------------------------------
//

#include <portinfo>

#include "IntegratorDomain.hh" // implementation of object methods

#include "pylith/feassemble/UpdateStateVars.hh" // HOLDSA UpdateStateVars
#include "pylith/feassemble/DSLabelAccess.hh" // USES DSLabelAccess
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
                                                         PetscFormKey key,
                                                         PetscIS cellIS,
                                                         PetscReal time,
                                                         PetscVec locX,
                                                         PetscVec locX_t,
                                                         PetscReal t,
                                                         PetscVec locF,
                                                         void *user);

extern "C" PetscErrorCode DMPlexComputeJacobian_Internal(PetscDM dm,
                                                         PetscFormKey key,
                                                         PetscIS cellIS,
                                                         PetscReal t,
                                                         PetscReal X_tShift,
                                                         PetscVec X,
                                                         PetscVec X_t,
                                                         PetscMat Jac,
                                                         PetscMat JacP,
                                                         void *user);

extern "C" PetscErrorCode DMPlexComputeJacobian_Action_Internal(PetscDM,
                                                                PetscFormKey,
                                                                PetscIS,
                                                                PetscReal,
                                                                PetscReal,
                                                                PetscVec,
                                                                PetscVec,
                                                                PetscVec,
                                                                PetscVec,
                                                                void *);

// ---------------------------------------------------------------------------------------------------------------------
// Default constructor.
pylith::feassemble::IntegratorDomain::IntegratorDomain(pylith::problems::Physics* const physics) :
    Integrator(physics),
    _materialMesh(NULL),
    _updateState(NULL) {
    GenericComponent::setName("integratordomain");
    _labelName = pylith::topology::Mesh::getCellsLabelName();
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
// Get mesh associated with integration domain.
const pylith::topology::Mesh&
pylith::feassemble::IntegratorDomain::getPhysicsDomainMesh(void) const {
    assert(_materialMesh);
    return *_materialMesh;
} // domainMesh


// ---------------------------------------------------------------------------------------------------------------------
void
pylith::feassemble::IntegratorDomain::setKernelsResidual(const std::vector<ResidualKernels>& kernels,
                                                         const pylith::topology::Field& solution) {
    PYLITH_METHOD_BEGIN;
    PYLITH_JOURNAL_DEBUG(_labelName<<"="<<_labelValue<<" setKernelsResidual(# kernels="<<kernels.size()<<")");

    PetscErrorCode err;
    DSLabelAccess dsLabel(solution.getDM(), _labelName.c_str(), _labelValue);
    for (size_t i = 0; i < kernels.size(); ++i) {
        const PetscInt i_field = solution.getSubfieldInfo(kernels[i].subfield.c_str()).index;
        const PetscInt i_part = kernels[i].part;
        err = PetscWeakFormAddResidual(dsLabel.weakForm(), dsLabel.label(), dsLabel.value(), i_field, i_part,
                                       kernels[i].r0, kernels[i].r1);PYLITH_CHECK_ERROR(err);

        switch (kernels[i].part) {
        case RESIDUAL_LHS:
            _hasLHSResidual = true;
            break;
        case RESIDUAL_RHS:
            _hasRHSResidual = true;
            break;
        default:
            PYLITH_JOURNAL_LOGICERROR("Unknown residual part " << kernels[i].part <<".");
        } // switch
    } // for

    pythia::journal::debug_t debug(GenericComponent::getName());
    if (debug.state()) {
        err = PetscDSView(dsLabel.ds(), PETSC_VIEWER_STDOUT_WORLD);PYLITH_CHECK_ERROR(err);
    } // if

    PYLITH_METHOD_END;
} // setKernelsResidual


// ---------------------------------------------------------------------------------------------------------------------
void
pylith::feassemble::IntegratorDomain::setKernelsJacobian(const std::vector<JacobianKernels>& kernels,
                                                         const pylith::topology::Field& solution) {
    PYLITH_METHOD_BEGIN;
    PYLITH_JOURNAL_DEBUG(_labelName<<"="<<_labelValue<<" setKernelsLHSJacobian(# kernels="<<kernels.size()<<")");

    PetscErrorCode err;
    DSLabelAccess dsLabel(solution.getDM(), _labelName.c_str(), _labelValue);
    for (size_t i = 0; i < kernels.size(); ++i) {
        const PetscInt i_fieldTrial = solution.getSubfieldInfo(kernels[i].subfieldTrial.c_str()).index;
        const PetscInt i_fieldBasis = solution.getSubfieldInfo(kernels[i].subfieldBasis.c_str()).index;
        const PetscInt i_part = kernels[i].part;
        err = PetscWeakFormAddJacobian(dsLabel.weakForm(), dsLabel.label(), dsLabel.value(), i_fieldTrial, i_fieldBasis,
                                       i_part, kernels[i].j0, kernels[i].j1, kernels[i].j2, kernels[i].j3);
        PYLITH_CHECK_ERROR(err);

        switch (kernels[i].part) {
        case JACOBIAN_LHS:
            _hasLHSJacobian = true;
            break;
        case JACOBIAN_LHS_LUMPED_INV:
            _hasLHSJacobianLumped = true;
            break;
        default:
            PYLITH_JOURNAL_LOGICERROR("Unknown Jacobian part " << kernels[i].part <<".");
        } // switch
    } // for

    pythia::journal::debug_t debug(GenericComponent::getName());
    if (debug.state()) {
        err = PetscDSView(dsLabel.ds(), PETSC_VIEWER_STDOUT_WORLD);PYLITH_CHECK_ERROR(err);
    } // if

    PYLITH_METHOD_END;
} // setKernelsLHSJacobian


// ---------------------------------------------------------------------------------------------------------------------
void
pylith::feassemble::IntegratorDomain::setKernelsUpdateStateVars(const std::vector<ProjectKernels>& kernels) {
    PYLITH_METHOD_BEGIN;
    PYLITH_JOURNAL_DEBUG(_labelName<<"="<<_labelValue<<" setKernelsUpdateStateVars(# kernels="<<kernels.size()<<")");

    _kernelsUpdateStateVars = kernels;

    PYLITH_METHOD_END;
} // setKernelsUpdateStateVars


// ---------------------------------------------------------------------------------------------------------------------
void
pylith::feassemble::IntegratorDomain::setKernelsDerivedField(const std::vector<ProjectKernels>& kernels) {
    PYLITH_METHOD_BEGIN;
    PYLITH_JOURNAL_DEBUG(_labelName<<"="<<_labelValue<<" setKernelsDerivedField(# kernels="<<kernels.size()<<")");

    _kernelsDerivedField = kernels;

    PYLITH_METHOD_END;
} // setKernelsDerivedField


// ---------------------------------------------------------------------------------------------------------------------
// Initialize integration domain, auxiliary field, and derived field. Update observers.
void
pylith::feassemble::IntegratorDomain::initialize(const pylith::topology::Field& solution) {
    PYLITH_METHOD_BEGIN;
    PYLITH_JOURNAL_DEBUG(_labelName<<"="<<_labelValue<<" intialize(solution="<<solution.getLabel()<<")");

    delete _materialMesh;
    _materialMesh = pylith::topology::MeshOps::createSubdomainMesh(solution.getMesh(), _labelName.c_str(), _labelValue, ":UNKOWN:");
    pylith::topology::CoordsVisitor::optimizeClosure(_materialMesh->getDM());

    Integrator::initialize(solution);

    assert(_auxiliaryField);
    PetscErrorCode err;
    PetscDM dmSoln = solution.getDM();assert(dmSoln);
    PetscDMLabel dmLabel = NULL;
    err = DMGetLabel(dmSoln, _labelName.c_str(), &dmLabel);PYLITH_CHECK_ERROR(err);
    err = DMSetAuxiliaryVec(dmSoln, dmLabel, _labelValue, _auxiliaryField->getLocalVector());PYLITH_CHECK_ERROR(err);

    if (_kernelsUpdateStateVars.size() > 0) {
        delete _updateState;_updateState = new pylith::feassemble::UpdateStateVars;assert(_updateState);
        _updateState->initialize(*_auxiliaryField);
    } // if

    pythia::journal::debug_t debug(GenericComponent::getName());
    if (debug.state()) {
        PYLITH_JOURNAL_DEBUG("Viewing auxiliary field.");
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
    PYLITH_JOURNAL_DEBUG(_labelName<<"="<<_labelValue<<" computeRHSResidual(residual="<<residual<<", t="<<t<<", dt="<<dt<<", solution="<<solution.getLabel()<<")");
    if (!_hasRHSResidual) { PYLITH_METHOD_END;}
    assert(residual);

    DSLabelAccess dsLabel(solution.getDM(), _labelName.c_str(), _labelValue);
    _setKernelConstants(solution, dt);

    PetscFormKey key;
    key.label = dsLabel.label();
    key.value = dsLabel.value();
    key.part = pylith::feassemble::Integrator::RESIDUAL_RHS;

    PetscErrorCode err;
    assert(solution.getLocalVector());
    assert(residual->getLocalVector());
    PetscVec solutionDotVec = NULL;
    err = DMPlexComputeResidual_Internal(dsLabel.dm(), key, dsLabel.cellsIS(), PETSC_MIN_REAL, solution.getLocalVector(),
                                         solutionDotVec, t, residual->getLocalVector(), NULL);PYLITH_CHECK_ERROR(err);

    PYLITH_METHOD_END;
} // computeRHSResidual


// ---------------------------------------------------------------------------------------------------------------------
// Compute LHS residual for F(t,s,\dot{s}).
void
pylith::feassemble::IntegratorDomain::computeLHSResidual(pylith::topology::Field* residual,
                                                         const PylithReal t,
                                                         const PylithReal dt,
                                                         const pylith::topology::Field& solution,
                                                         const pylith::topology::Field& solutionDot) {
    PYLITH_METHOD_BEGIN;
    PYLITH_JOURNAL_DEBUG(_labelName<<"="<<_labelValue<<" computeLHSResidual(residual="<<residual<<", t="<<t<<", dt="<<dt<<", solution="<<solution.getLabel()<<")");
    if (!_hasLHSResidual) { PYLITH_METHOD_END; }

    DSLabelAccess dsLabel(solution.getDM(), _labelName.c_str(), _labelValue);
    _setKernelConstants(solution, dt);

    PetscFormKey key;
    key.label = dsLabel.label();
    key.value = dsLabel.value();
    key.part = pylith::feassemble::Integrator::RESIDUAL_LHS;

    PetscErrorCode err;
    assert(solution.getLocalVector());
    assert(solutionDot.getLocalVector());
    assert(residual->getLocalVector());
    err = DMPlexComputeResidual_Internal(dsLabel.dm(), key, dsLabel.cellsIS(), PETSC_MIN_REAL, solution.getLocalVector(),
                                         solutionDot.getLocalVector(), t, residual->getLocalVector(), NULL);PYLITH_CHECK_ERROR(err);

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
    PYLITH_JOURNAL_DEBUG(_labelName<<"="<<_labelValue<<" computeLHSJacobian(jacobianMat="<<jacobianMat<<", precondMat="<<precondMat<<", t="<<t<<", dt="<<dt<<", solution="<<solution.getLabel()<<", solutionDot="<<solutionDot.getLabel()<<")");

    _needNewLHSJacobian = false;
    if (!_hasLHSJacobian) { PYLITH_METHOD_END;}

    DSLabelAccess dsLabel(solution.getDM(), _labelName.c_str(), _labelValue);
    _setKernelConstants(solution, dt);

    PetscFormKey key;
    key.label = dsLabel.label();
    key.value = dsLabel.value();
    key.part = pylith::feassemble::Integrator::JACOBIAN_LHS;

    PetscErrorCode err;
    assert(solution.getLocalVector());
    assert(solutionDot.getLocalVector());
    assert(jacobianMat);
    assert(precondMat);
    err = DMPlexComputeJacobian_Internal(dsLabel.dm(), key, dsLabel.cellsIS(), t, s_tshift, solution.getLocalVector(),
                                         solutionDot.getLocalVector(), jacobianMat, precondMat, NULL);PYLITH_CHECK_ERROR(err);

    PYLITH_METHOD_END;
} // computeLHSJacobian


// ---------------------------------------------------------------------------------------------------------------------
// Compute inverse of lumped LHS Jacobian for F(t,s,\dot{s}).
void
pylith::feassemble::IntegratorDomain::computeLHSJacobianLumpedInv(pylith::topology::Field* jacobianInv,
                                                                  const PylithReal t,
                                                                  const PylithReal dt,
                                                                  const PylithReal s_tshift,
                                                                  const pylith::topology::Field& solution) {
    PYLITH_METHOD_BEGIN;
    PYLITH_JOURNAL_DEBUG(_labelName<<"="<<_labelValue<<" computeLHSJacobianLumpedInv(jacobianInv="<<jacobianInv<<", t="<<t<<", dt="<<dt<<", solution="<<solution.getLabel()<<")");

    _needNewLHSJacobianLumped = false;
    if (!_hasLHSJacobianLumped) { PYLITH_METHOD_END;}

    DSLabelAccess dsLabel(solution.getDM(), _labelName.c_str(), _labelValue);
    _setKernelConstants(solution, dt);

    PetscFormKey key;
    key.label = dsLabel.label();
    key.value = dsLabel.value();
    key.part = pylith::feassemble::Integrator::JACOBIAN_LHS_LUMPED_INV;

    PetscErrorCode err;
    PetscVec vecRowSum = NULL;
    err = DMGetLocalVector(dsLabel.dm(), &vecRowSum);PYLITH_CHECK_ERROR(err);
    err = VecSet(vecRowSum, 1.0);PYLITH_CHECK_ERROR(err);

    assert(jacobianInv);
    assert(jacobianInv->getLocalVector());
    err = DMPlexComputeJacobian_Action_Internal(dsLabel.dm(), key, dsLabel.cellsIS(), t, s_tshift, vecRowSum, NULL,
                                                vecRowSum, jacobianInv->getLocalVector(), NULL);PYLITH_CHECK_ERROR(err);

    // Compute the Jacobian inverse.
    err = VecReciprocal(jacobianInv->getLocalVector());PYLITH_CHECK_ERROR(err);

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
    PetscDMLabel dmLabel = NULL;
    PetscInt labelValue = 0;
    err = DMSetAuxiliaryVec(stateVarsDM, dmLabel, labelValue, _auxiliaryField->getLocalVector());PYLITH_CHECK_ERROR(err);
    err = DMProjectFieldLocal(stateVarsDM, t, solution.getLocalVector(), kernelsStateVars, INSERT_VALUES,
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
        const pylith::topology::Field::SubfieldInfo& sinfo = _derivedField->getSubfieldInfo(_kernelsDerivedField[iKernel].subfield.c_str());
        kernelsArray[sinfo.index] = _kernelsDerivedField[iKernel].f;
    } // for

    PetscErrorCode err = 0;

    PetscDM derivedDM = _derivedField->getDM();
    assert(_auxiliaryField);
    PetscDMLabel dmLabel = NULL;
    PetscInt labelValue = 0;
    err = DMSetAuxiliaryVec(derivedDM, dmLabel, labelValue, _auxiliaryField->getLocalVector());PYLITH_CHECK_ERROR(err);
    err = DMProjectFieldLocal(derivedDM, t, solution.getLocalVector(), kernelsArray, INSERT_VALUES, _derivedField->getLocalVector());PYLITH_CHECK_ERROR(err);
    delete[] kernelsArray;kernelsArray = NULL;

    pythia::journal::debug_t debug(GenericComponent::getName());
    if (debug.state()) {
        PYLITH_JOURNAL_DEBUG("Viewing derived field.");
        _derivedField->view("Derived field");
    } // if

    PYLITH_METHOD_END;
} // _computeDerivedField


// End of file
