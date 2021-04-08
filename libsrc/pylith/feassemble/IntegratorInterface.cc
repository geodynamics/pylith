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

#include "pylith/feassemble/InterfacePatches.hh" // USES InterfacePatches
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

            /** Set cohesive cell integration kernels.
             *
             * @param[in] integrator Integrator for interface.
             * @param[in] kernels Integration kernels (pointwise functions) for cohesive cells.
             * @param[in] solution Field with current trial solution.
             */
            static
            void setWeakFormKernels(const pylith::feassemble::IntegratorInterface* integrator,
                                    const std::vector<pylith::feassemble::IntegratorInterface::ResidualKernels>& kernels,
                                    const pylith::topology::Field& solution);

            /** Set cohesive cell integration kernels.
             *
             * @param[in] integrator Integrator for interface.
             * @param[in] kernels Integration kernels (pointwise functions) for cohesive cells.
             * @param[in] solution Field with current trial solution.
             */
            static
            void transferWeakFormKernels(const pylith::feassemble::IntegratorInterface* integrator,
                                         const pylith::topology::Field& solution);

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
    _interfaceSurfaceLabel(""),
    _integrationPatches(NULL) {
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
pylith::feassemble::IntegratorInterface::setIntegrationPatches(pylith::feassemble::InterfacePatches* patches) {
    PYLITH_METHOD_BEGIN;
    PYLITH_JOURNAL_DEBUG("setIntegrationPatches(patches="<<patches<<")");

    delete _integrationPatches;_integrationPatches = patches;

    PYLITH_METHOD_END;
} // setIntegrationPatches


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
    typedef InterfacePatches::keysmap_t keysmap_t;

    pythia::journal::debug_t debug(_IntegratorInterface::genericComponent);
    debug << pythia::journal::at(__HERE__)
          << "_IntegratorInterface::computeResidual(residual="<<typeid(residual).name()<<", integrator"<<typeid(integrator).name()
          <<"# kernels="<<kernels.size()<<", t="<<t<<", dt="<<dt<<", solution="<<solution.getLabel()<<", solutionDot="
          <<solutionDot.getLabel()<<")"
          << pythia::journal::endl;

    assert(integrator);
    assert(residual);

    setWeakFormKernels(integrator, kernels, solution);
    transferWeakFormKernels(integrator, solution);

    // Loop over integration patches.
    PetscErrorCode err = 0;
    PetscDM dmSoln = solution.dmMesh();
    const pylith::topology::Field* auxiliaryField = integrator->getAuxiliaryField();assert(auxiliaryField);
    const InterfacePatches* patches = integrator->_integrationPatches;assert(patches);
    const keysmap_t& keysmap = patches->getKeys();
    for (keysmap_t::const_iterator iter = keysmap.begin(); iter != keysmap.end(); ++iter) {
        PetscHashFormKey weakFormKeys[3];
        weakFormKeys[0] = iter->second.negative.petscKey(solution);
        weakFormKeys[1] = iter->second.positive.petscKey(solution);
        weakFormKeys[2] = iter->second.cohesive.petscKey(solution);

        const PetscInt labelValue = weakFormKeys[2].value;
        PetscIS cohesiveCellIS = NULL;
        PetscInt numCohesiveCells = 0;
        const PetscInt* cohesiveCells = NULL;
        err = DMGetStratumIS(dmSoln, integrator->getLabelName(), labelValue, &cohesiveCellIS);PYLITH_CHECK_ERROR(err);
        err = ISGetSize(cohesiveCellIS, &numCohesiveCells);PYLITH_CHECK_ERROR(err);assert(numCohesiveCells > 0);
        err = ISGetIndices(cohesiveCellIS, &cohesiveCells);PYLITH_CHECK_ERROR(err);assert(cohesiveCells);
        assert(pylith::topology::MeshOps::isCohesiveCell(dmSoln, cohesiveCells[0]));

        // Get auxiliary data
        // :TODO: FIX THIS. This needs to be updated. We need to provide the auxiliary field(s) for the cohesive cells
        // and the adjacent cells on the negative and positive sides of the fault.
        PetscErrorCode err = DMSetAuxiliaryVec(dmSoln, NULL, 0, auxiliaryField->localVector());PYLITH_CHECK_ERROR(err);

        assert(solution.localVector());
        assert(residual->localVector());
        err = DMPlexComputeResidual_Hybrid_Internal(dmSoln, weakFormKeys, cohesiveCellIS, t, solution.localVector(),
                                                    solutionDot.localVector(), t,
                                                    residual->localVector(), NULL);PYLITH_CHECK_ERROR(err);
        err = ISRestoreIndices(cohesiveCellIS, &cohesiveCells);PYLITH_CHECK_ERROR(err);
        err = ISDestroy(&cohesiveCellIS);PYLITH_CHECK_ERROR(err);
    } // for

#if 0
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
#endif

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
pylith::feassemble::_IntegratorInterface::setWeakFormKernels(const pylith::feassemble::IntegratorInterface* integrator,
                                                             const std::vector<pylith::feassemble::IntegratorInterface::ResidualKernels>& kernels,
                                                             const pylith::topology::Field& solution) {
    PYLITH_METHOD_BEGIN;
    typedef InterfacePatches::keysmap_t keysmap_t;

    assert(integrator);
    const InterfacePatches* patches = integrator->_integrationPatches;assert(patches);

    PetscErrorCode err = 0;
    const keysmap_t& keysmap = patches->getKeys();
    for (keysmap_t::const_iterator iter = keysmap.begin(); iter != keysmap.end(); ++iter) {
        PetscHashFormKey key = iter->second.cohesive.petscKey(solution);
        const PetscWeakForm weakForm = iter->second.cohesive.getWeakForm();

        for (size_t i = 0; i < kernels.size(); ++i) {
            err = PetscWeakFormSetIndexBdResidual(weakForm, key.label, key.value, key.field,
                                                  i, kernels[i].r0, i, kernels[i].r1);PYLITH_CHECK_ERROR(err);
        } // for
    } // for

    { // DEBUGGING
        PetscIS cohesiveCells = NULL;
        PetscInt numCohesiveCells = 0;
        const PetscInt* cellIndices = NULL;
        err = DMGetStratumIS(solution.dmMesh(), integrator->getLabelName(), integrator->getLabelValue(), &cohesiveCells);PYLITH_CHECK_ERROR(err);
        err = ISGetSize(cohesiveCells, &numCohesiveCells);PYLITH_CHECK_ERROR(err);assert(numCohesiveCells > 0);
        err = ISGetIndices(cohesiveCells, &cellIndices);PYLITH_CHECK_ERROR(err);assert(cellIndices);

        assert(pylith::topology::MeshOps::isCohesiveCell(solution.dmMesh(), cellIndices[0]));
        PetscDS prob = NULL;
        err = DMGetCellDS(solution.dmMesh(), cellIndices[0], &prob);PYLITH_CHECK_ERROR(err);
        err = PetscDSView(prob, PETSC_VIEWER_STDOUT_WORLD);
    } // DEBUGGING

    PYLITH_METHOD_END;
} // TMP_setWeakFormKernels


// ---------------------------------------------------------------------------------------------------------------------
void
pylith::feassemble::_IntegratorInterface::transferWeakFormKernels(const pylith::feassemble::IntegratorInterface* integrator,
                                                                  const pylith::topology::Field& solution) {
    PYLITH_METHOD_BEGIN;
    typedef InterfacePatches::keysmap_t keysmap_t;

    assert(integrator);
    const InterfacePatches* patches = integrator->_integrationPatches;assert(patches);

    PetscErrorCode err = 0;
    const keysmap_t& keysmap = patches->getKeys();
    for (keysmap_t::const_iterator iter = keysmap.begin(); iter != keysmap.end(); ++iter) {
        PetscBdPointFunc* f0 = NULL;
        PetscBdPointFunc* f1 = NULL;
        PylithInt f0Count = 0, f1Count = 0;

        const PetscWeakForm weakFormCohesive = iter->second.cohesive.getWeakForm();

        const PetscWeakForm weakFormNegative = iter->second.negative.getWeakForm();
        if (weakFormNegative) {
            PetscHashFormKey keyNegative = iter->second.negative.petscKey(solution);
            err = PetscWeakFormGetBdResidual(weakFormNegative, keyNegative.label, keyNegative.value, keyNegative.field, &f0Count, &f0, &f1Count, &f1);PYLITH_CHECK_ERROR(err);
            err = PetscWeakFormSetBdResidual(weakFormCohesive, keyNegative.label, keyNegative.value, keyNegative.field, f0Count, f0, f1Count, f1);PYLITH_CHECK_ERROR(err);
        } // if

        const PetscWeakForm weakFormPositive = iter->second.positive.getWeakForm();
        if (weakFormPositive) {
            PetscHashFormKey keyPositive = iter->second.positive.petscKey(solution);
            err = PetscWeakFormGetBdResidual(weakFormPositive, keyPositive.label, keyPositive.value, keyPositive.field, &f0Count, &f0, &f1Count, &f1);PYLITH_CHECK_ERROR(err);
            err = PetscWeakFormSetBdResidual(weakFormCohesive, keyPositive.label, keyPositive.value, keyPositive.field, f0Count, f0, f1Count, f1);PYLITH_CHECK_ERROR(err);
        } // if
    } // for

    PYLITH_METHOD_END;
} // transferWeakForms


// End of file
