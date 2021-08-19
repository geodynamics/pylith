// -*- C++ -*-
//
// ------------------------------------------------------------------------------------------------
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
// ------------------------------------------------------------------------------------------------
//

#include <portinfo>

#include "IntegratorInterface.hh" // implementation of object methods

#include "pylith/feassemble/InterfacePatches.hh" // USES InterfacePatches
#include "pylith/feassemble/DSLabelAccess.hh" // USES DSLabelAccess
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

// ------------------------------------------------------------------------------------------------
// Local "private" functions.
namespace pylith {
    namespace feassemble {
        class _IntegratorInterface {
public:

            /** Compute residual using current kernels.
             *
             * @param[out] residual Field for residual.
             * @param[in] integrator Integrator for boundary.
             * @param[in] residualPart Residual part to compute.
             * @param[in] t Current time.
             * @param[in] dt Current time step.
             * @param[in] solution Field with current trial solution.
             * @param[in] solutionDot Field with time derivative of current trial solution.
             */
            static
            void computeResidual(pylith::topology::Field* residual,
                                 const pylith::feassemble::IntegratorInterface* integrator,
                                 pylith::feassemble::Integrator::ResidualPart residualPart,
                                 const PylithReal t,
                                 const PylithReal dt,
                                 const pylith::topology::Field& solution,
                                 const pylith::topology::Field& solutionDot);

            /** Compute Jacobian using current kernels.
             *
             * @param[out] jacobianMat PETSc Mat with Jacobian sparse matrix.
             * @param[out] precondMat PETSc Mat with Jacobian preconditioning sparse matrix.
             * @param[in] integrator Integrator for boundary.
             * @param[in] jacobianPart Jacobian part to compute.
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
                                 pylith::feassemble::Integrator::JacobianPart jacobianPart,
                                 const PylithReal t,
                                 const PylithReal dt,
                                 const PylithReal s_tshift,
                                 const pylith::topology::Field& solution,
                                 const pylith::topology::Field& solutionDot);

            static const char* genericComponent;

        }; // _IntegratorInterface
        const char* _IntegratorInterface::genericComponent = "integratorinterface";

    } // feassemble
} // pylith

// ------------------------------------------------------------------------------------------------
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


// ------------------------------------------------------------------------------------------------
// Destructor.
pylith::feassemble::IntegratorInterface::~IntegratorInterface(void) {
    deallocate();
} // destructor


// ------------------------------------------------------------------------------------------------
// Deallocate PETSc and local data structures.
void
pylith::feassemble::IntegratorInterface::deallocate(void) {
    PYLITH_METHOD_BEGIN;

    pylith::feassemble::Integrator::deallocate();

    delete _interfaceMesh;_interfaceMesh = NULL;
    delete _integrationPatches;_integrationPatches = NULL;

    PYLITH_METHOD_END;
} // deallocate


// ------------------------------------------------------------------------------------------------
// Set label marking boundary associated with boundary condition surface.
void
pylith::feassemble::IntegratorInterface::setSurfaceMarkerLabel(const char* value) {
    PYLITH_JOURNAL_DEBUG("setSurfaceMarkerLabel(value="<<value<<")");

    if (strlen(value) == 0) {
        throw std::runtime_error("Empty string given for boundary condition label.");
    } // if

    _interfaceSurfaceLabel = value;
} // setSurfaceMarkerLabel


// ------------------------------------------------------------------------------------------------
// Get label marking boundary associated with boundary condition surface.
const char*
pylith::feassemble::IntegratorInterface::getSurfaceMarkerLabel(void) const {
    return _interfaceSurfaceLabel.c_str();
} // getSurfaceMarkerLabel


// ------------------------------------------------------------------------------------------------
// Get mesh associated with integration domain.
const pylith::topology::Mesh&
pylith::feassemble::IntegratorInterface::getPhysicsDomainMesh(void) const {
    assert(_interfaceMesh);
    return *_interfaceMesh;
} // domainMesh


// ------------------------------------------------------------------------------------------------
// Set weak form keys for integration patch.
void
pylith::feassemble::IntegratorInterface::setIntegrationPatches(pylith::feassemble::InterfacePatches* patches) {
    PYLITH_METHOD_BEGIN;
    PYLITH_JOURNAL_DEBUG("setIntegrationPatches(patches="<<patches<<")");

    delete _integrationPatches;_integrationPatches = patches;

    PYLITH_METHOD_END;
} // setIntegrationPatches


// ------------------------------------------------------------------------------------------------
// Set kernels for residual.
void
pylith::feassemble::IntegratorInterface::setKernelsResidual(const std::vector<ResidualKernels>& kernels,
                                                            const pylith::topology::Field& solution) {
    PYLITH_METHOD_BEGIN;
    PYLITH_JOURNAL_DEBUG(_labelName<<"="<<_labelValue<<" setKernelsResidual(# kernels="<<kernels.size()<<")");
    typedef InterfacePatches::keysmap_t keysmap_t;

    PetscErrorCode err;
    const keysmap_t& keysmap = _integrationPatches->getKeys();
    for (keysmap_t::const_iterator iter = keysmap.begin(); iter != keysmap.end(); ++iter) {
        const PetscWeakForm weakForm = iter->second.cohesive.getWeakForm();

        for (size_t i = 0; i < kernels.size(); ++i) {
            PetscFormKey key;
            switch (kernels[i].face) {
            case IntegratorInterface::NEGATIVE_FACE:
                key = iter->second.negative.getPetscKey(solution, kernels[i].part, kernels[i].subfield.c_str());
                break;
            case IntegratorInterface::POSITIVE_FACE:
                key = iter->second.positive.getPetscKey(solution, kernels[i].part, kernels[i].subfield.c_str());
                break;
            case IntegratorInterface::FAULT_FACE:
                key = iter->second.cohesive.getPetscKey(solution, kernels[i].part, kernels[i].subfield.c_str());
                break;
            default:
                PYLITH_JOURNAL_LOGICERROR("Unknown integration face.");
            } // switch
            err = PetscWeakFormAddBdResidual(weakForm, key.label, key.value, key.field, key.part,
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
    } // for

    pythia::journal::debug_t debug(_IntegratorInterface::genericComponent);
    if (debug.state()) {
        DSLabelAccess dsLabel(solution.getDM(), _labelName.c_str(), _labelValue);
        err = PetscDSView(dsLabel.ds(), PETSC_VIEWER_STDOUT_WORLD);PYLITH_CHECK_ERROR(err);
    } // if

    PYLITH_METHOD_END;
} // setKernelsRHSResidual


// ------------------------------------------------------------------------------------------------
// Set kernels for Jacobian.
void
pylith::feassemble::IntegratorInterface::setKernelsJacobian(const std::vector<JacobianKernels>& kernels,
                                                            const pylith::topology::Field& solution) {
    PYLITH_METHOD_BEGIN;
    PYLITH_JOURNAL_DEBUG(_labelName<<"="<<_labelValue<<" setKernelsJacobian(# kernels="<<kernels.size()<<")");
    typedef InterfacePatches::keysmap_t keysmap_t;

    PetscErrorCode err = 0;
    const keysmap_t& keysmap = _integrationPatches->getKeys();
    for (keysmap_t::const_iterator iter = keysmap.begin(); iter != keysmap.end(); ++iter) {
        const PetscWeakForm weakForm = iter->second.cohesive.getWeakForm();

        PetscInt numFields = 0;
        err = PetscWeakFormGetNumFields(weakForm, &numFields);PYLITH_CHECK_ERROR(err);

        for (size_t i = 0; i < kernels.size(); ++i) {
            PetscFormKey key;
            switch (kernels[i].face) {
            case IntegratorInterface::NEGATIVE_FACE:
                key = iter->second.negative.getPetscKey(solution, kernels[i].part, kernels[i].subfieldTrial.c_str(),
                                                        kernels[i].subfieldBasis.c_str());
                break;
            case IntegratorInterface::POSITIVE_FACE:
                key = iter->second.positive.getPetscKey(solution, kernels[i].part, kernels[i].subfieldTrial.c_str(),
                                                        kernels[i].subfieldBasis.c_str());
                break;
            case IntegratorInterface::FAULT_FACE:
                key = iter->second.cohesive.getPetscKey(solution, kernels[i].part, kernels[i].subfieldTrial.c_str(),
                                                        kernels[i].subfieldBasis.c_str());
                break;
            default:
                PYLITH_JOURNAL_LOGICERROR("Unknown integration face.");
            } // switch
            const PetscInt i_trial = key.field / numFields;
            const PetscInt i_basis = key.field % numFields;
            err = PetscWeakFormAddBdJacobian(weakForm, key.label, key.value, i_trial, i_basis, key.part,
                                             kernels[i].j0, kernels[i].j1, kernels[i].j2, kernels[i].j3);PYLITH_CHECK_ERROR(err);

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
    } // for

    pythia::journal::debug_t debug(_IntegratorInterface::genericComponent);
    if (debug.state()) {
        DSLabelAccess dsLabel(solution.getDM(), _labelName.c_str(), _labelValue);
        err = PetscDSView(dsLabel.ds(), PETSC_VIEWER_STDOUT_WORLD);PYLITH_CHECK_ERROR(err);
    } // if

    PYLITH_METHOD_END;
} // setKernelsJacobian


// ------------------------------------------------------------------------------------------------
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

    PetscErrorCode err = 0;
    PetscDM dmSoln = solution.getDM();
    typedef InterfacePatches::keysmap_t keysmap_t;
    const keysmap_t& keysmap = _integrationPatches->getKeys();
    const pylith::topology::Field* auxiliaryField = getAuxiliaryField();assert(auxiliaryField);
    for (keysmap_t::const_iterator iter = keysmap.begin(); iter != keysmap.end(); ++iter) {
        PetscFormKey key = iter->second.cohesive.getPetscKey(solution, pylith::feassemble::Integrator::RESIDUAL_LHS);
        err = DMSetAuxiliaryVec(dmSoln, key.label, key.value, auxiliaryField->getLocalVector());PYLITH_CHECK_ERROR(err);
    } // for

    PYLITH_METHOD_END;
} // initialize


// ------------------------------------------------------------------------------------------------
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


// ------------------------------------------------------------------------------------------------
// Compute RHS residual for G(t,s).
void
pylith::feassemble::IntegratorInterface::computeRHSResidual(pylith::topology::Field* residual,
                                                            const PylithReal t,
                                                            const PylithReal dt,
                                                            const pylith::topology::Field& solution) {
    PYLITH_METHOD_BEGIN;
    PYLITH_JOURNAL_DEBUG("computeRHSResidual(residual="<<residual<<", t="<<t<<", dt="<<dt<<", solution="<<solution.getLabel()<<")");
    if (!_hasRHSResidual) { PYLITH_METHOD_END;}

    _setKernelConstants(solution, dt);

    pylith::topology::Field solutionDot(solution.getMesh()); // No dependence on time derivative of solution in RHS.
    solutionDot.setLabel("solution_dot");

    _IntegratorInterface::computeResidual(residual, this, pylith::feassemble::Integrator::RESIDUAL_RHS,
                                          t, dt, solution, solutionDot);

    PYLITH_METHOD_END;
} // computeRHSResidual


// ------------------------------------------------------------------------------------------------
// Compute LHS residual for F(t,s,\dot{s}).
void
pylith::feassemble::IntegratorInterface::computeLHSResidual(pylith::topology::Field* residual,
                                                            const PylithReal t,
                                                            const PylithReal dt,
                                                            const pylith::topology::Field& solution,
                                                            const pylith::topology::Field& solutionDot) {
    PYLITH_METHOD_BEGIN;
    PYLITH_JOURNAL_DEBUG("computeLHSResidual(residual="<<residual<<", t="<<t<<", dt="<<dt<<", solution="<<solution.getLabel()<<")");
    if (!_hasLHSResidual) { PYLITH_METHOD_END;}

    _setKernelConstants(solution, dt);

    _IntegratorInterface::computeResidual(residual, this, pylith::feassemble::Integrator::RESIDUAL_LHS,
                                          t, dt, solution, solutionDot);

    PYLITH_METHOD_END;
} // computeLHSResidual


// ------------------------------------------------------------------------------------------------
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
    if (!_hasLHSJacobian) { PYLITH_METHOD_END;}

    _setKernelConstants(solution, dt);

    _IntegratorInterface::computeJacobian(jacobianMat, precondMat, this, pylith::feassemble::Integrator::JACOBIAN_LHS,
                                          t, dt, s_tshift, solution, solutionDot);

    PYLITH_METHOD_END;
} // computeLHSJacobian


// ------------------------------------------------------------------------------------------------
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


// ------------------------------------------------------------------------------------------------
// Compute residual.
void
pylith::feassemble::_IntegratorInterface::computeResidual(pylith::topology::Field* residual,
                                                          const pylith::feassemble::IntegratorInterface* integrator,
                                                          pylith::feassemble::Integrator::ResidualPart residualPart,
                                                          const PylithReal t,
                                                          const PylithReal dt,
                                                          const pylith::topology::Field& solution,
                                                          const pylith::topology::Field& solutionDot) {
    PYLITH_METHOD_BEGIN;
    typedef InterfacePatches::keysmap_t keysmap_t;

    pythia::journal::debug_t debug(_IntegratorInterface::genericComponent);
    debug << pythia::journal::at(__HERE__)
          << "_IntegratorInterface::computeResidual(residual="<<typeid(residual).name()<<", integrator"<<typeid(integrator).name()
          <<", residualPart="<<residualPart<<", t="<<t<<", dt="<<dt<<", solution="<<solution.getLabel()<<", solutionDot="
          <<solutionDot.getLabel()<<")"
          << pythia::journal::endl;

    assert(integrator);
    assert(residual);

    // Loop over integration patches.
    PetscErrorCode err = 0;
    PetscDM dmSoln = solution.getDM();
    const InterfacePatches* patches = integrator->_integrationPatches;assert(patches);
    const keysmap_t& keysmap = patches->getKeys();
    for (keysmap_t::const_iterator iter = keysmap.begin(); iter != keysmap.end(); ++iter) {
        PetscFormKey weakFormKeys[3];
        weakFormKeys[0] = iter->second.negative.getPetscKey(solution, residualPart, "displacement");
        weakFormKeys[1] = iter->second.positive.getPetscKey(solution, residualPart, "displacement");
        weakFormKeys[2] = iter->second.cohesive.getPetscKey(solution, residualPart, "lagrange_multiplier_fault");

        PetscIS patchCellsIS = NULL;
        PetscInt numPatchCells = 0;
        const PetscInt* patchCells = NULL;
        err = DMGetStratumIS(dmSoln, patches->getLabelName(), weakFormKeys[2].value, &patchCellsIS);PYLITH_CHECK_ERROR(err);
        err = ISGetSize(patchCellsIS, &numPatchCells);PYLITH_CHECK_ERROR(err);assert(numPatchCells > 0);
        err = ISGetIndices(patchCellsIS, &patchCells);PYLITH_CHECK_ERROR(err);assert(patchCells);
        assert(pylith::topology::MeshOps::isCohesiveCell(dmSoln, patchCells[0]));

        assert(solution.getLocalVector());
        assert(residual->getLocalVector());
        err = DMPlexComputeResidual_Hybrid_Internal(dmSoln, weakFormKeys, patchCellsIS, t, solution.getLocalVector(),
                                                    solutionDot.getLocalVector(), t,
                                                    residual->getLocalVector(), NULL);PYLITH_CHECK_ERROR(err);
        err = ISRestoreIndices(patchCellsIS, &patchCells);PYLITH_CHECK_ERROR(err);
        err = ISDestroy(&patchCellsIS);PYLITH_CHECK_ERROR(err);
    } // for

    PYLITH_METHOD_END;
} // computeResidual


// ------------------------------------------------------------------------------------------------
// Compute Jacobian.
void
pylith::feassemble::_IntegratorInterface::computeJacobian(PetscMat jacobianMat,
                                                          PetscMat precondMat,
                                                          const pylith::feassemble::IntegratorInterface* integrator,
                                                          pylith::feassemble::Integrator::JacobianPart jacobianPart,
                                                          const PylithReal t,
                                                          const PylithReal dt,
                                                          const PylithReal s_tshift,
                                                          const pylith::topology::Field& solution,
                                                          const pylith::topology::Field& solutionDot) {
    PYLITH_METHOD_BEGIN;
    typedef InterfacePatches::keysmap_t keysmap_t;

    pythia::journal::debug_t debug(_IntegratorInterface::genericComponent);
    debug << pythia::journal::at(__HERE__)
          << "_IntegratorInterface::computeJacobian(jacobianMat="<<jacobianMat<<", precondMat"<<precondMat
          <<", integrator"<<typeid(integrator).name()<<", jacobianPart="<<jacobianPart<<", t="<<t<<", dt="<<dt<<", solution="
          <<solution.getLabel()<<", solutionDot="<<solutionDot.getLabel()<<")"
          << pythia::journal::endl;

    assert(jacobianMat);
    assert(precondMat);

    PetscErrorCode err;
    PetscDM dmSoln = solution.getDM();
    const InterfacePatches* patches = integrator->_integrationPatches;assert(patches);
    const keysmap_t& keysmap = patches->getKeys();
    for (keysmap_t::const_iterator iter = keysmap.begin(); iter != keysmap.end(); ++iter) {
        PetscFormKey weakFormKeys[3];
        weakFormKeys[0] = iter->second.negative.getPetscKey(solution, jacobianPart, "displacement", "lagrange_multiplier_fault");
        weakFormKeys[1] = iter->second.positive.getPetscKey(solution, jacobianPart, "displacement", "lagrange_multiplier_fault");
        weakFormKeys[2] = iter->second.cohesive.getPetscKey(solution, jacobianPart, "lagrange_multiplier_fault", "displacement");

        PetscIS patchCellsIS = NULL;
        PetscInt numPatchCells = 0;
        const PetscInt* patchCells = NULL;
        err = DMGetStratumIS(dmSoln, patches->getLabelName(), weakFormKeys[2].value, &patchCellsIS);PYLITH_CHECK_ERROR(err);
        err = ISGetSize(patchCellsIS, &numPatchCells);PYLITH_CHECK_ERROR(err);assert(numPatchCells > 0);
        err = ISGetIndices(patchCellsIS, &patchCells);PYLITH_CHECK_ERROR(err);assert(patchCells);
        assert(pylith::topology::MeshOps::isCohesiveCell(dmSoln, patchCells[0]));

        assert(solution.getLocalVector());
        err = DMPlexComputeJacobian_Hybrid_Internal(dmSoln, weakFormKeys, patchCellsIS, t, s_tshift, solution.getLocalVector(),
                                                    solutionDot.getLocalVector(), jacobianMat, precondMat,
                                                    NULL);PYLITH_CHECK_ERROR(err);
        err = ISRestoreIndices(patchCellsIS, &patchCells);PYLITH_CHECK_ERROR(err);
        err = ISDestroy(&patchCellsIS);PYLITH_CHECK_ERROR(err);
    }
    PYLITH_METHOD_END;
} // computeJacobian


// End of file
