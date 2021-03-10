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

#include "IntegratorInterface.hh" // implementation of object methods

#include "pylith/problems/Physics.hh" // USES Physics
#include "pylith/topology/Mesh.hh" // USES Mesh
#include "pylith/topology/Field.hh" // USES Field
#include "pylith/topology/CoordsVisitor.hh" // USES CoordsVisitor::optimizeClosure()
#include "pylith/faults/TopologyOps.hh" // USES TopologyOps
#include "pylith/topology/MeshOps.hh" // USES MeshOps

#include "spatialdata/spatialdb/GravityField.hh" // HASA GravityField
#include "petscds.h" // USES PetscDS

#include "pylith/utils/journals.hh" // USES PYLITH_JOURNAL_*
#include <cassert> // USES assert()
#include <typeinfo> // USES typeid()
#include <stdexcept> // USES std::runtime_error

extern "C" PetscErrorCode DMPlexComputeResidual_Hybrid_Internal(PetscDM dm,
                                                                PetscFormKey key[],
                                                                PetscIS cellIS,
                                                                PetscReal time,
                                                                PetscVec locX,
                                                                PetscVec locX_t,
                                                                PetscReal t,
                                                                PetscVec locF,
                                                                void *user);

extern "C" PetscErrorCode DMPlexComputeJacobian_Hybrid_Internal(PetscDM dm,
                                                                PetscFormKey key[],
                                                                PetscIS cellIS,
                                                                PetscReal t,
                                                                PetscReal X_tShift,
                                                                PetscVec locX,
                                                                PetscVec locX_t,
                                                                PetscMat Jac,
                                                                PetscMat JacP,
                                                                void *user);

// ---------------------------------------------------------------------------------------------------------------------
// Local "private" functions.
namespace pylith {
    namespace feassemble {
        class _IntegratorInterface {
public:

            /** Compute residual using current kernels.
             *
             * @param[out] residual Field for residual.
             * @param[in] integrator Integrator for boundary.
             * @param[in] kernels Kernels for computing residual.
             * @param[in] t Current time.
             * @param[in] dt Current time step.
             * @param[in] solution Field with current trial solution.
             * @param[in] solutionDot Field with time derivative of current trial solution.
             */
            static
            void computeResidual(pylith::topology::Field* residual,
                                 const pylith::feassemble::IntegratorInterface* integrator,
                                 const std::vector<pylith::feassemble::IntegratorInterface::ResidualKernels>& kernels,
                                 const PylithReal t,
                                 const PylithReal dt,
                                 const pylith::topology::Field& solution,
                                 const pylith::topology::Field& solutionDot);

            /** Compute Jacobian using current kernels.
             *
             * @param[out] jacobianMat PETSc Mat with Jacobian sparse matrix.
             * @param[out] precondMat PETSc Mat with Jacobian preconditioning sparse matrix.
             * @param[in] integrator Integrator for boundary.
             * @param[in] kernels Kernels for computing Jacobian.
             * @param[in] t Current time.
             * @param[in] dt Current time step.
             * @param[in] s_tshift Scale for time derivative.
             * @param[in] solution Field with current trial solution.
             * @param[in] solutionDot Field with time derivative of current trial solution.
             */
            static
            void computeJacobian(PetscMat jacobianMat,
                                 PetscMat precondMat,
                                 const pylith::feassemble::IntegratorInterface* integrator,
                                 const std::vector<pylith::feassemble::IntegratorInterface::JacobianKernels>& kernels,
                                 const PylithReal t,
                                 const PylithReal dt,
                                 const PylithReal s_tshift,
                                 const pylith::topology::Field& solution,
                                 const pylith::topology::Field& solutionDot);

            static
            void transferWeakForms(const pylith::topology::Field& solution,
                                   const pylith::feassemble::IntegratorInterface::WeakFormKeys& weakFormKeys,
                                   const PetscInt cohesiveCell,
                                   const PetscInt adjacentCellNegative,
                                   const PetscInt adjacentCellPositive);

            static const char* genericComponent;
        }; // _IntegratorInterface
        const char* _IntegratorInterface::genericComponent = "integratorinterface";

    } // feassemble
} // pylith

// ---------------------------------------------------------------------------------------------------------------------
// Default constructor.
pylith::feassemble::IntegratorInterface::IntegratorInterface(pylith::problems::Physics* const physics) :
    Integrator(physics),
    _interfaceMesh(NULL),
    _interfaceSurfaceLabel("") {
    GenericComponent::setName(_IntegratorInterface::genericComponent);
    _labelValue = 100;
    _labelName = pylith::topology::Mesh::getCellsLabelName();
} // constructor


// ---------------------------------------------------------------------------------------------------------------------
// Destructor.
pylith::feassemble::IntegratorInterface::~IntegratorInterface(void) {
    deallocate();
} // destructor


// ---------------------------------------------------------------------------------------------------------------------
// Deallocate PETSc and local data structures.
void
pylith::feassemble::IntegratorInterface::deallocate(void) {
    PYLITH_METHOD_BEGIN;

    pylith::feassemble::Integrator::deallocate();

    delete _interfaceMesh;_interfaceMesh = NULL;

    PYLITH_METHOD_END;
} // deallocate


// ---------------------------------------------------------------------------------------------------------------------
// Set label marking boundary associated with boundary condition surface.
void
pylith::feassemble::IntegratorInterface::setSurfaceMarkerLabel(const char* value) {
    PYLITH_JOURNAL_DEBUG("setSurfaceMarkerLabel(value="<<value<<")");

    if (strlen(value) == 0) {
        throw std::runtime_error("Empty string given for boundary condition label.");
    } // if

    _interfaceSurfaceLabel = value;
} // setSurfaceMarkerLabel


// ---------------------------------------------------------------------------------------------------------------------
// Get label marking boundary associated with boundary condition surface.
const char*
pylith::feassemble::IntegratorInterface::getSurfaceMarkerLabel(void) const {
    return _interfaceSurfaceLabel.c_str();
} // getSurfaceMarkerLabel


// ---------------------------------------------------------------------------------------------------------------------
// Get mesh associated with integration domain.
const pylith::topology::Mesh&
pylith::feassemble::IntegratorInterface::getPhysicsDomainMesh(void) const {
    assert(_interfaceMesh);
    return *_interfaceMesh;
} // domainMesh


// ---------------------------------------------------------------------------------------------------------------------
// Set weak form keys for integration patch.
void
pylith::feassemble::IntegratorInterface::setPatchWeakFormKeys(const int labelValue,
                                                              const WeakFormKeys& keys) {
    PYLITH_METHOD_BEGIN;
    PYLITH_JOURNAL_DEBUG("setPatchWeakFormKeys(labelValue="<<labelValue<<", keys="<<typeid(keys).name()<<")");

    PYLITH_JOURNAL_LOGICERROR(":TODO: Implement.");

    PYLITH_METHOD_END;
} // setPatchWeakFormKeys


// ---------------------------------------------------------------------------------------------------------------------
// Set kernels for RHS residual.
void
pylith::feassemble::IntegratorInterface::setKernelsRHSResidual(const std::vector<ResidualKernels>& kernels) {
    PYLITH_METHOD_BEGIN;
    PYLITH_JOURNAL_DEBUG("setKernelsRHSResidual(# kernels="<<kernels.size()<<")");

    _kernelsRHSResidual = kernels;

    PYLITH_METHOD_END;
} // setKernelsRHSResidual


// ---------------------------------------------------------------------------------------------------------------------
// Set kernels for LHS residual.
void
pylith::feassemble::IntegratorInterface::setKernelsLHSResidual(const std::vector<ResidualKernels>& kernels) {
    PYLITH_METHOD_BEGIN;
    PYLITH_JOURNAL_DEBUG("setKernelsLHSResidual(# kernels="<<kernels.size()<<")");

    _kernelsLHSResidual = kernels;

    PYLITH_METHOD_END;
} // setKernelsLHSResidual


// ---------------------------------------------------------------------------------------------------------------------
// Set kernels for LHS Jacobian.
void
pylith::feassemble::IntegratorInterface::setKernelsLHSJacobian(const std::vector<JacobianKernels>& kernels) {
    PYLITH_METHOD_BEGIN;
    PYLITH_JOURNAL_DEBUG("setKernelsLHSJacobian(# kernels="<<kernels.size()<<")");

    _kernelsLHSJacobian = kernels;

    PYLITH_METHOD_END;
} // setKernelsLHSJacobian


// ---------------------------------------------------------------------------------------------------------------------
// Initialize integration domain, auxiliary field, and derived field. Update observers.
void
pylith::feassemble::IntegratorInterface::initialize(const pylith::topology::Field& solution) {
    PYLITH_METHOD_BEGIN;
    PYLITH_JOURNAL_DEBUG("intialize(solution="<<solution.getLabel()<<")");

    const bool isSubmesh = true;
    delete _interfaceMesh;_interfaceMesh = new pylith::topology::Mesh(isSubmesh);assert(_interfaceMesh);
    pylith::faults::TopologyOps::createFaultParallel(_interfaceMesh, solution.getMesh(), _labelValue, _labelName.c_str(),
                                                     _interfaceSurfaceLabel.c_str());
    pylith::topology::MeshOps::checkTopology(*_interfaceMesh);
    pylith::topology::CoordsVisitor::optimizeClosure(_interfaceMesh->getDM());

    Integrator::initialize(solution);

    PYLITH_METHOD_END;
} // initialize


// ---------------------------------------------------------------------------------------------------------------------
// Update auxiliary field values to current time.
void
pylith::feassemble::IntegratorInterface::updateState(const double t) {
    PYLITH_METHOD_BEGIN;
    PYLITH_JOURNAL_DEBUG("updateState(t="<<t<<")");

    Integrator::updateState(t);

    assert(_physics);
    _physics->updateAuxiliaryField(_auxiliaryField, t);

    pythia::journal::debug_t debug(GenericComponent::getName());
    if (debug.state()) {
        assert(_auxiliaryField);
        PYLITH_JOURNAL_DEBUG("IntegratorInterface component '" << GenericComponent::getName() << "' for '"
                                                               <<_physics->getIdentifier()
                                                               << "': viewing auxiliary field.");
        _auxiliaryField->view("IntegratorInterface auxiliary field", pylith::topology::Field::VIEW_ALL);
    } // if

    PYLITH_METHOD_END;
} // updateState


// ---------------------------------------------------------------------------------------------------------------------
// Compute RHS residual for G(t,s).
void
pylith::feassemble::IntegratorInterface::computeRHSResidual(pylith::topology::Field* residual,
                                                            const PylithReal t,
                                                            const PylithReal dt,
                                                            const pylith::topology::Field& solution) {
    PYLITH_METHOD_BEGIN;
    PYLITH_JOURNAL_DEBUG("computeRHSResidual(residual="<<residual<<", t="<<t<<", dt="<<dt<<", solution="<<solution.getLabel()<<")");

    if (0 == _kernelsRHSResidual.size()) { PYLITH_METHOD_END;}

    _setKernelConstants(solution, dt);

    pylith::topology::Field solutionDot(solution.getMesh()); // No dependence on time derivative of solution in RHS.
    solutionDot.setLabel("solution_dot");

    _IntegratorInterface::computeResidual(residual, this, _kernelsRHSResidual, t, dt, solution, solutionDot);

    PYLITH_METHOD_END;
} // computeRHSResidual


// ---------------------------------------------------------------------------------------------------------------------
// Compute LHS residual for F(t,s,\dot{s}).
void
pylith::feassemble::IntegratorInterface::computeLHSResidual(pylith::topology::Field* residual,
                                                            const PylithReal t,
                                                            const PylithReal dt,
                                                            const pylith::topology::Field& solution,
                                                            const pylith::topology::Field& solutionDot) {
    PYLITH_METHOD_BEGIN;
    PYLITH_JOURNAL_DEBUG("computeLHSResidual(residual="<<residual<<", t="<<t<<", dt="<<dt<<", solution="<<solution.getLabel()<<")");

    if (0 == _kernelsLHSResidual.size()) { PYLITH_METHOD_END;}

    _setKernelConstants(solution, dt);

    _IntegratorInterface::computeResidual(residual, this, _kernelsLHSResidual, t, dt, solution, solutionDot);

    PYLITH_METHOD_END;
} // computeLHSResidual


// ---------------------------------------------------------------------------------------------------------------------
// Compute LHS Jacobian for F(t,s,\dot{s}).
void
pylith::feassemble::IntegratorInterface::computeLHSJacobian(PetscMat jacobianMat,
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

    _IntegratorInterface::computeJacobian(jacobianMat, precondMat, this, _kernelsLHSJacobian, t, dt, s_tshift,
                                          solution, solutionDot);

    PYLITH_METHOD_END;
} // computeLHSJacobian


// ---------------------------------------------------------------------------------------------------------------------
// Compute LHS Jacobian for F(t,s,\dot{s}).
void
pylith::feassemble::IntegratorInterface::computeLHSJacobianLumpedInv(pylith::topology::Field* jacobianInv,
                                                                     const PylithReal t,
                                                                     const PylithReal dt,
                                                                     const PylithReal s_tshift,
                                                                     const pylith::topology::Field& solution) {
    PYLITH_METHOD_BEGIN;
    PYLITH_JOURNAL_DEBUG("computeLHSJacobianLumpedInv(jacobianInv="<<jacobianInv<<", t="<<t<<", dt="<<dt<<", solution="<<solution.getLabel()<<") empty method");

    _needNewLHSJacobianLumped = false;
    // No implementation needed for interface.

    PYLITH_METHOD_END;
} // computeLHSJacobianLumpedInv


// ---------------------------------------------------------------------------------------------------------------------
// Compute residual.
void
pylith::feassemble::_IntegratorInterface::computeResidual(pylith::topology::Field* residual,
                                                          const pylith::feassemble::IntegratorInterface* integrator,
                                                          const std::vector<pylith::feassemble::IntegratorInterface::ResidualKernels>& kernels,
                                                          const PylithReal t,
                                                          const PylithReal dt,
                                                          const pylith::topology::Field& solution,
                                                          const pylith::topology::Field& solutionDot) {
    PYLITH_METHOD_BEGIN;
    pythia::journal::debug_t debug(_IntegratorInterface::genericComponent);
    debug << pythia::journal::at(__HERE__)
          << "_IntegratorInterface::computeResidual(residual="<<typeid(residual).name()<<", integrator"<<typeid(integrator).name()
          <<"# kernels="<<kernels.size()<<", t="<<t<<", dt="<<dt<<", solution="<<solution.getLabel()<<", solutionDot="
          <<solutionDot.getLabel()<<")"
          << pythia::journal::endl;

    assert(integrator);
    assert(residual);

    const pylith::topology::Field* auxiliaryField = integrator->getAuxiliaryField();assert(auxiliaryField);

    // Get auxiliary data
    PetscDMLabel dmLabel = NULL;
    PetscInt labelValue = 0;
    PetscDM dmSoln = solution.dmMesh();
    PetscErrorCode err = DMSetAuxiliaryVec(dmSoln, dmLabel, labelValue, auxiliaryField->localVector());PYLITH_CHECK_ERROR(err);

    // Loop over integration patches using label.
    PetscIS labelValueIS = NULL;
    PetscInt numLabelValues = 0;
    const PetscInt* labelValues = NULL;
    err = DMGetLabelIdIS(dmSoln, integrator->getLabelName(), &labelValueIS);PYLITH_CHECK_ERROR(err);
    err = ISGetSize(labelValueIS, &numLabelValues);PYLITH_CHECK_ERROR(err);assert(numLabelValues > 0);
    err = ISGetIndices(labelValueIS, &labelValues);PYLITH_CHECK_ERROR(err);assert(labelValues);
    for (PetscInt iValue = 0; iValue < numLabelValues; ++iValue) {
        const PetscInt labelValue = labelValues[iValue];

        PetscIS cohesiveCellIS = NULL;
        PetscInt numCohesiveCells = 0;
        const PetscInt* cohesiveCells = NULL;
        err = DMGetStratumIS(dmSoln, integrator->getLabelName(), labelValue, &cohesiveCellIS);PYLITH_CHECK_ERROR(err);
        err = ISGetSize(cohesiveCellIS, &numCohesiveCells);PYLITH_CHECK_ERROR(err);assert(numCohesiveCells > 0);
        err = ISGetIndices(cohesiveCellIS, &cohesiveCells);PYLITH_CHECK_ERROR(err);assert(cohesiveCells);
        assert(pylith::topology::MeshOps::isCohesiveCell(dmSoln, cohesiveCells[0]));

        PetscFormKey weakFormKeys[3];
        weakFormKeys[0].label = NULL;
        weakFormKeys[0].value = 0;
        weakFormKeys[0].field = 0;
        weakFormKeys[0].part = 0;
        weakFormKeys[1].label = NULL;
        weakFormKeys[1].value = 0;
        weakFormKeys[1].field = 0;
        weakFormKeys[1].part = 0;
        weakFormKeys[2].label = NULL;
        weakFormKeys[2].value = 0;
        weakFormKeys[2].field = 0;
        weakFormKeys[2].part = 0;

        // Get auxiliary data
        err = DMSetAuxiliaryVec(dmSoln, NULL, 0, auxiliaryField->localVector());PYLITH_CHECK_ERROR(err);

        assert(solution.localVector());
        assert(residual->localVector());
        err = DMPlexComputeResidual_Hybrid_Internal(dmSoln, weakFormKeys, cohesiveCellIS, t, solution.localVector(),
                                                    solutionDot.localVector(), t,
                                                    residual->localVector(), NULL);PYLITH_CHECK_ERROR(err);
        err = ISRestoreIndices(cohesiveCellIS, &cohesiveCells);PYLITH_CHECK_ERROR(err);
        err = ISDestroy(&cohesiveCellIS);PYLITH_CHECK_ERROR(err);
    } // for
    err = ISRestoreIndices(labelValueIS, &labelValues);PYLITH_CHECK_ERROR(err);
    err = ISDestroy(&labelValueIS);PYLITH_CHECK_ERROR(err);

#if 0
    PetscDS prob = NULL;
    err = DMGetCellDS(dmSoln, cohesiveCells[0], &prob);PYLITH_CHECK_ERROR(err);
    // Compute the local residual
    for (size_t i = 0; i < kernels.size(); ++i) {
        const PetscInt i_field = solution.getSubfieldInfo(kernels[i].subfield.c_str()).index;
        err = PetscDSSetBdResidual(prob, i_field, kernels[i].r0, kernels[i].r1);PYLITH_CHECK_ERROR(err);
    } // for
    err = DMPlexComputeResidual_Hybrid_Internal(dmSoln, keys, cohesiveCellsIS, t, solution.getLocalVector(),
                                                solutionDot.getLocalVector(), t,
                                                residual->getLocalVector(), NULL);PYLITH_CHECK_ERROR(err);
    err = ISRestoreIndices(cohesiveCellsIS, &cellIndices);PYLITH_CHECK_ERROR(err);
    err = ISDestroy(&cohesiveCellsIS);PYLITH_CHECK_ERROR(err);
#endif

    PYLITH_METHOD_END;
} // computeResidual


// ---------------------------------------------------------------------------------------------------------------------
// Compute Jacobian.
void
pylith::feassemble::_IntegratorInterface::computeJacobian(PetscMat jacobianMat,
                                                          PetscMat precondMat,
                                                          const pylith::feassemble::IntegratorInterface* integrator,
                                                          const std::vector<pylith::feassemble::IntegratorInterface::JacobianKernels>& kernels,
                                                          const PylithReal t,
                                                          const PylithReal dt,
                                                          const PylithReal s_tshift,
                                                          const pylith::topology::Field& solution,
                                                          const pylith::topology::Field& solutionDot) {
    PYLITH_METHOD_BEGIN;
    pythia::journal::debug_t debug(_IntegratorInterface::genericComponent);
    debug << pythia::journal::at(__HERE__)
          << "_IntegratorInterface::computeJacobian(jacobianMat="<<jacobianMat<<", precondMat"<<precondMat
          <<", integrator"<<typeid(integrator).name()<<"# kernels="<<kernels.size()<<", t="<<t<<", dt="<<dt<<", solution="
          <<solution.getLabel()<<", solutionDot="<<solutionDot.getLabel()<<")"
          << pythia::journal::endl;

    assert(jacobianMat);
    assert(precondMat);

    const pylith::topology::Field* auxiliaryField = integrator->getAuxiliaryField();assert(auxiliaryField);

    PetscErrorCode err;

    PetscDM dmSoln = solution.getDM();

    PetscIS cohesiveCells = NULL;
    PetscInt numCohesiveCells = 0;
    const PetscInt* cellIndices = NULL;
    err = DMGetStratumIS(dmSoln, integrator->getLabelName(), integrator->getLabelValue(), &cohesiveCells);PYLITH_CHECK_ERROR(err);
    err = ISGetSize(cohesiveCells, &numCohesiveCells);PYLITH_CHECK_ERROR(err);assert(numCohesiveCells > 0);
    err = ISGetIndices(cohesiveCells, &cellIndices);PYLITH_CHECK_ERROR(err);assert(cellIndices);

    assert(pylith::topology::MeshOps::isCohesiveCell(dmSoln, cellIndices[0]));
    PetscDS prob = NULL;
    err = DMGetCellDS(dmSoln, cellIndices[0], &prob);PYLITH_CHECK_ERROR(err);

    // Set auxiliary data
    PetscDMLabel dmLabel = NULL;
    PetscInt labelValue = 0;
    err = DMSetAuxiliaryVec(dmSoln, dmLabel, labelValue, auxiliaryField->getLocalVector());PYLITH_CHECK_ERROR(err);

    PetscFormKey keys[3];
    keys[0].label = NULL;
    keys[0].value = 0;
    keys[0].field = 0;
    keys[0].part = 0;
    keys[1].label = NULL;
    keys[1].value = 0;
    keys[1].field = 0;
    keys[1].part = 0;
    keys[2].label = NULL;
    keys[2].value = 0;
    keys[2].field = 0;
    keys[2].part = 0;

    // Compute the local Jacobian
    for (size_t i = 0; i < kernels.size(); ++i) {
        const PetscInt i_fieldTrial = solution.getSubfieldInfo(kernels[i].subfieldTrial.c_str()).index;
        const PetscInt i_fieldBasis = solution.getSubfieldInfo(kernels[i].subfieldBasis.c_str()).index;
        err = PetscDSSetBdJacobian(prob, i_fieldTrial, i_fieldBasis, kernels[i].j0, kernels[i].j1, kernels[i].j2, kernels[i].j3);PYLITH_CHECK_ERROR(err);
    } // for

    assert(solution.getLocalVector());
    err = DMPlexComputeJacobian_Hybrid_Internal(dmSoln, keys, cohesiveCells, t, s_tshift, solution.getLocalVector(),
                                                solutionDot.getLocalVector(), jacobianMat, precondMat,
                                                NULL);PYLITH_CHECK_ERROR(err);
    err = ISRestoreIndices(cohesiveCells, &cellIndices);PYLITH_CHECK_ERROR(err);
    err = ISDestroy(&cohesiveCells);PYLITH_CHECK_ERROR(err);

    PYLITH_METHOD_END;
} // computeJacobian


// ---------------------------------------------------------------------------------------------------------------------
void
pylith::feassemble::_IntegratorInterface::transferWeakForms(const pylith::topology::Field& solution,
                                                            const pylith::feassemble::IntegratorInterface::WeakFormKeys& weakFormKeys,
                                                            const PetscInt cohesiveCell,
                                                            const PetscInt adjacentCellNegative,
                                                            const PetscInt adjacentCellPositive) {
    PYLITH_METHOD_BEGIN;
#if 0

    const char* const patchLabelName = weakFormKeys.cohesive.name.c_str();

    const char* const materialsLabelName = weakFormKeys.negative.name.c_str();
    assert(weakFormKeys.negative.name == weakFormKeys.positive.name);
    const PylithInt labelValueNegative = weakFormKeys.negative.value;
    const PylithInt labelValuePositive = weakFormKeys.positive.value;

    PetscDM dmSoln = solution.dmMesh();
    const PetscInt i_lagrange = solution.subfieldInfo("lagrange_multiplier_fault").index;

    // Transfer residual kernels for negative/positive face from material weak form to fault weak form.
    PetscErrorCode err = 0;
    PetscDMLabel labelMaterial = NULL, labelCohesive = NULL;
    err = DMGetLabel(dmSoln, materialsLabelName, &labelMaterial);PYLITH_CHECK_ERROR(err);
    err = DMGetLabel(dmSoln, patchLabelName, &labelCohesive);PYLITH_CHECK_ERROR(err);

    PetscDS probMaterial = NULL, probCohesive = NULL;
    PetscWeakForm weakFormMaterial = NULL, weakFormCohesive = NULL;
    err = DMGetCellDS(dmSoln, cohesiveCell, &probCohesive);PYLITH_CHECK_ERROR(err);
    err = PetscDSGetWeakForm(probCohesive, &weakFormCohesive);PYLITH_CHECK_ERROR(err);

    PetscBdPointFunc* f0 = NULL;
    PetscBdPointFunc* f1 = NULL;
    PylithInt f0Count = 0, f1Count = 0;
    err = DMGetCellDS(dmSoln, adjacentCellNegative, &probMaterial);PYLITH_CHECK_ERROR(err);
    err = PetscDSGetWeakForm(probMaterial, &weakFormMaterial);PYLITH_CHECK_ERROR(err);
    err = PetscWeakFormGetBdResidual(weakFormMaterial, labelMaterial, labelValueNegative, i_lagrange, &f0Count, &f0, &f1Count, &f1);PYLITH_CHECK_ERROR(err);
    err = PetscWeakFormSetBdResidual(weakFormCohesive, labelMaterial, labelValueNegative, i_lagrange, f0Count, f0, f1Count, f1);PYLITH_CHECK_ERROR(err);

    err = DMGetCellDS(dmSoln, adjacentCellPositive, &probMaterial);PYLITH_CHECK_ERROR(err);
    err = PetscDSGetWeakForm(probMaterial, &weakFormMaterial);PYLITH_CHECK_ERROR(err);
    err = PetscWeakFormGetBdResidual(weakFormMaterial, labelMaterial, labelValuePositive, i_lagrange, &f0Count, &f0, &f1Count, &f1);PYLITH_CHECK_ERROR(err);
    err = PetscWeakFormSetBdResidual(weakFormCohesive, labelCohesive, labelValuePositive, i_lagrange, f0Count, f0, f1Count, f1);PYLITH_CHECK_ERROR(err);
#endif

    PYLITH_METHOD_END;
} // _transferWeakForms


// End of file
