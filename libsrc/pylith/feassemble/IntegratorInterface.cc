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
#include "pylith/problems/IntegrationData.hh" // USES IntegrationData
#include "pylith/topology/Mesh.hh" // USES Mesh
#include "pylith/topology/Field.hh" // USES Field
#include "pylith/topology/CoordsVisitor.hh" // USES CoordsVisitor::optimizeClosure()
#include "pylith/faults/TopologyOps.hh" // USES TopologyOps
#include "pylith/topology/MeshOps.hh" // USES MeshOps
#include "pylith/materials/Material.hh" // USES Material

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
            typedef pylith::feassemble::IntegratorInterface::ResidualKernels ResidualKernels;
            typedef pylith::feassemble::IntegratorInterface::JacobianKernels JacobianKernels;
public:

            /** Compute weak form key part for face.
             *
             * For integration with hybrid cells, we must distinguish among integration on the
             * negative face, positive face, and fault face as well as the residual term.
             */
            static
            PetscInt computeWeakFormPart(const PetscInt part,
                                         const PetscInt face);

            /** Compute residual using current kernels.
             *
             * @param[out] residual Field for residual.
             * @param[in] integrator Integrator for boundary.
             * @param[in] residualPart Residual part to compute.
             * @param[in] integrationData Data needed to integrate governing equations.
             */
            static
            void computeResidual(pylith::topology::Field* residual,
                                 const pylith::feassemble::IntegratorInterface* integrator,
                                 pylith::feassemble::Integrator::ResidualPart residualPart,
                                 const pylith::problems::IntegrationData& integrationData);

            /** Compute Jacobian using current kernels.
             *
             * @param[out] jacobianMat PETSc Mat with Jacobian sparse matrix.
             * @param[out] precondMat PETSc Mat with Jacobian preconditioning sparse matrix.
             * @param[in] integrator Integrator for boundary.
             * @param[in] jacobianPart Jacobian part to compute.
             * @param[in] integrationData Data needed to integrate governing equations.
             */
            static
            void computeJacobian(PetscMat jacobianMat,
                                 PetscMat precondMat,
                                 const pylith::feassemble::IntegratorInterface* integrator,
                                 pylith::feassemble::Integrator::JacobianPart jacobianPart,
                                 const pylith::problems::IntegrationData& integrationData);

            /** Add bounding material kernels to array of residual kernels.
             *
             * @param[inout] kernels Array of residual kernels.
             * @param[in] face Side of interior interface for integration.
             * @param[in] key Key associated with integration patch for bounding material.
             * @param[in] solution Field with current trial solution.
             * @param[in] materials Materials in model.
             */
            static
            void addMaterialKernels(std::vector<ResidualKernels>* kernels,
                                    const pylith::feassemble::IntegratorInterface::FaceEnum face,
                                    const pylith::feassemble::FEKernelKey key,
                                    const pylith::topology::Field& solution,
                                    const std::vector<pylith::materials::Material*>& materials);

            /** Add bounding material kernels to array of Jacobian kernels.
             *
             * @param[inout] kernels Array of Jacobian kernels.
             * @param[in] face Side of interior interface for integration.
             * @param[in] key Key associated with integration patch for bounding material.
             * @param[in] solution Field with current trial solution.
             * @param[in] materials Materials in model.
             */
            static
            void addMaterialKernels(std::vector<JacobianKernels>* kernels,
                                    const pylith::feassemble::IntegratorInterface::FaceEnum face,
                                    const pylith::feassemble::FEKernelKey key,
                                    const pylith::topology::Field& solution,
                                    const std::vector<pylith::materials::Material*>& materials);

            /** Get material corresponding to label name and value.
             *
             * @param[in] materials Materials in model.
             * @param[in] labelName Name of label.
             * @param[in] labelValue Label value.
             * @returns Material with label name and value.
             */
            static
            pylith::materials::Material* getMaterial(const std::vector<pylith::materials::Material*> materials,
                                                     const char* labelName,
                                                     const int labelValue);

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
pylith::feassemble::IntegratorInterface::setKernels(const std::vector<ResidualKernels>& kernels,
                                                    const pylith::topology::Field& solution,
                                                    const std::vector<pylith::materials::Material*>& materials) {
    PYLITH_METHOD_BEGIN;
    PYLITH_JOURNAL_DEBUG(_labelName<<"="<<_labelValue<<" setKernels(Residual)(# kernels="<<kernels.size()<<")");
    typedef InterfacePatches::keysmap_t keysmap_t;

    PetscErrorCode err;
    const keysmap_t& keysmap = _integrationPatches->getKeys();
    for (keysmap_t::const_iterator iter = keysmap.begin(); iter != keysmap.end(); ++iter) {
        const PetscWeakForm weakForm = iter->second.cohesive.getWeakForm();
        std::vector<ResidualKernels> kernelsPatch(kernels);

#if 0
        _IntegratorInterface::addMaterialKernels(&kernelsPatch, pylith::feassemble::IntegratorInterface::NEGATIVE_FACE,
                                                 iter->second.negative, solution, materials);
        _IntegratorInterface::addMaterialKernels(&kernelsPatch, pylith::feassemble::IntegratorInterface::POSITIVE_FACE,
                                                 iter->second.positive, solution, materials);
#endif
        for (size_t i = 0; i < kernelsPatch.size(); ++i) {
            PetscFormKey key;

            const PetscInt interfacePart =
                _IntegratorInterface::computeWeakFormPart(kernelsPatch[i].part, kernelsPatch[i].face);

            switch (kernelsPatch[i].face) {
            case IntegratorInterface::NEGATIVE_FACE:
                key = iter->second.negative.getPetscKey(solution, interfacePart, kernelsPatch[i].subfield.c_str());
                break;
            case IntegratorInterface::POSITIVE_FACE:
                key = iter->second.positive.getPetscKey(solution, interfacePart, kernelsPatch[i].subfield.c_str());
                break;
            case IntegratorInterface::COHESIVE_FACE:
                key = iter->second.cohesive.getPetscKey(solution, interfacePart, kernelsPatch[i].subfield.c_str());
                break;
            default:
                PYLITH_JOURNAL_LOGICERROR("Unknown integration face ("<<kernelsPatch[i].face<<").");
            } // switch
            err = PetscWeakFormAddBdResidual(weakForm, key.label, key.value, key.field, key.part,
                                             kernelsPatch[i].r0, kernelsPatch[i].r1);PYLITH_CHECK_ERROR(err);

            switch (kernelsPatch[i].part) {
            case RESIDUAL_LHS:
                _hasLHSResidual = true;
                break;
            case RESIDUAL_RHS:
                _hasRHSResidual = true;
                break;
            default:
                PYLITH_JOURNAL_LOGICERROR("Unknown residual part " << kernelsPatch[i].part <<".");
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
pylith::feassemble::IntegratorInterface::setKernels(const std::vector<JacobianKernels>& kernels,
                                                    const pylith::topology::Field& solution,
                                                    const std::vector<pylith::materials::Material*>& materials) {
    PYLITH_METHOD_BEGIN;
    PYLITH_JOURNAL_DEBUG(_labelName<<"="<<_labelValue<<" setKernels(Jacobian)(# kernels="<<kernels.size()<<")");
    typedef InterfacePatches::keysmap_t keysmap_t;

    PetscErrorCode err = 0;
    const keysmap_t& keysmap = _integrationPatches->getKeys();
    for (keysmap_t::const_iterator iter = keysmap.begin(); iter != keysmap.end(); ++iter) {
        const PetscWeakForm weakForm = iter->second.cohesive.getWeakForm();
        std::vector<JacobianKernels> kernelsPatch(kernels);

        _IntegratorInterface::addMaterialKernels(&kernelsPatch, pylith::feassemble::IntegratorInterface::NEGATIVE_FACE,
                                                 iter->second.negative, solution, materials);
        _IntegratorInterface::addMaterialKernels(&kernelsPatch, pylith::feassemble::IntegratorInterface::POSITIVE_FACE,
                                                 iter->second.positive, solution, materials);

        PetscInt numFields = 0;
        err = PetscWeakFormGetNumFields(weakForm, &numFields);PYLITH_CHECK_ERROR(err);

        for (size_t i = 0; i < kernelsPatch.size(); ++i) {
            PetscFormKey key;
            const PetscInt interfacePart =
                _IntegratorInterface::computeWeakFormPart(kernelsPatch[i].part, kernelsPatch[i].face);

            switch (kernelsPatch[i].face) {
            case IntegratorInterface::NEGATIVE_FACE:
                key = iter->second.negative.getPetscKey(solution, interfacePart, kernelsPatch[i].subfieldTrial.c_str(),
                                                        kernelsPatch[i].subfieldBasis.c_str());
                break;
            case IntegratorInterface::POSITIVE_FACE:
                key = iter->second.positive.getPetscKey(solution, interfacePart, kernelsPatch[i].subfieldTrial.c_str(),
                                                        kernelsPatch[i].subfieldBasis.c_str());
                break;
            case IntegratorInterface::COHESIVE_FACE:
                key = iter->second.cohesive.getPetscKey(solution, interfacePart, kernelsPatch[i].subfieldTrial.c_str(),
                                                        kernelsPatch[i].subfieldBasis.c_str());
                break;
            default:
                PYLITH_JOURNAL_LOGICERROR("Unknown integration face (" << kernelsPatch[i].face << ").");
            } // switch
            const PetscInt i_trial = key.field / numFields;
            const PetscInt i_basis = key.field % numFields;
            err = PetscWeakFormAddBdJacobian(weakForm, key.label, key.value, i_trial, i_basis, key.part,
                                             kernelsPatch[i].j0, kernelsPatch[i].j1, kernelsPatch[i].j2, kernelsPatch[i].j3);PYLITH_CHECK_ERROR(err);

            switch (kernelsPatch[i].part) {
            case JACOBIAN_LHS:
                _hasLHSJacobian = true;
                break;
            case JACOBIAN_LHS_LUMPED_INV:
                _hasLHSJacobianLumped = true;
                break;
            default:
                PYLITH_JOURNAL_LOGICERROR("Unknown Jacobian part (" << kernelsPatch[i].part <<").");
            } // switch
        } // for
    } // for

    pythia::journal::debug_t debug(_IntegratorInterface::genericComponent);
    if (debug.state()) {
        DSLabelAccess dsLabel(solution.getDM(), _labelName.c_str(), _labelValue);
        err = PetscDSView(dsLabel.ds(), PETSC_VIEWER_STDOUT_WORLD);PYLITH_CHECK_ERROR(err);
    } // if

    PYLITH_METHOD_END;
} // setKernels


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
// Set auxiliary field values for current time.
void
pylith::feassemble::IntegratorInterface::setState(const PylithReal t,
                                                  const PylithReal dt) {
    PYLITH_METHOD_BEGIN;
    PYLITH_JOURNAL_DEBUG("setState(t="<<t<<", dt="<<dt<<")");

    Integrator::setState(t, dt);

    assert(_physics);
    _physics->updateAuxiliaryField(_auxiliaryField, t, dt);

    pythia::journal::debug_t debug(GenericComponent::getName());
    if (debug.state()) {
        assert(_auxiliaryField);
        PYLITH_JOURNAL_DEBUG("IntegratorInterface component '" << GenericComponent::getName() << "' for '"
                                                               <<_physics->getIdentifier()
                                                               << "': viewing auxiliary field.");
        _auxiliaryField->view("IntegratorInterface auxiliary field", pylith::topology::Field::VIEW_ALL);
    } // if

    PYLITH_METHOD_END;
} // setState


// ------------------------------------------------------------------------------------------------
// Compute RHS residual for G(t,s).
void
pylith::feassemble::IntegratorInterface::computeRHSResidual(pylith::topology::Field* residual,
                                                            const pylith::problems::IntegrationData& integrationData) {
    PYLITH_METHOD_BEGIN;
    PYLITH_JOURNAL_DEBUG("computeRHSResidual(residual="<<residual<<", integrationData="<<integrationData.str()<<")");
    if (!_hasRHSResidual) { PYLITH_METHOD_END;}

    _IntegratorInterface::computeResidual(residual, this, pylith::feassemble::Integrator::RESIDUAL_RHS, integrationData);

    PYLITH_METHOD_END;
} // computeRHSResidual


// ------------------------------------------------------------------------------------------------
// Compute LHS residual for F(t,s,\dot{s}).
void
pylith::feassemble::IntegratorInterface::computeLHSResidual(pylith::topology::Field* residual,
                                                            const pylith::problems::IntegrationData& integrationData) {
    PYLITH_METHOD_BEGIN;
    PYLITH_JOURNAL_DEBUG("computeLHSResidual(residual="<<residual<<", integrationData="<<integrationData.str()<<")");
    if (!_hasLHSResidual) { PYLITH_METHOD_END;}

    _IntegratorInterface::computeResidual(residual, this, pylith::feassemble::Integrator::RESIDUAL_LHS, integrationData);

    PYLITH_METHOD_END;
} // computeLHSResidual


// ------------------------------------------------------------------------------------------------
// Compute LHS Jacobian for F(t,s,\dot{s}).
void
pylith::feassemble::IntegratorInterface::computeLHSJacobian(PetscMat jacobianMat,
                                                            PetscMat precondMat,
                                                            const pylith::problems::IntegrationData& integrationData) {
    PYLITH_METHOD_BEGIN;
    PYLITH_JOURNAL_DEBUG("computeLHSJacobian(jacobianMat="<<jacobianMat<<", precondMat="<<precondMat<<", integrationData="<<integrationData.str()<<")");

    _needNewLHSJacobian = false;
    if (!_hasLHSJacobian) { PYLITH_METHOD_END;}

    _IntegratorInterface::computeJacobian(jacobianMat, precondMat, this, pylith::feassemble::Integrator::JACOBIAN_LHS, integrationData);

    PYLITH_METHOD_END;
} // computeLHSJacobian


// ------------------------------------------------------------------------------------------------
// Compute LHS Jacobian for F(t,s,\dot{s}).
void
pylith::feassemble::IntegratorInterface::computeLHSJacobianLumpedInv(pylith::topology::Field* jacobianInv,
                                                                     const pylith::problems::IntegrationData& integrationData) {
    PYLITH_METHOD_BEGIN;
    PYLITH_JOURNAL_DEBUG("computeLHSJacobianLumpedInv(jacobianInv="<<jacobianInv<<", integrationData="<<integrationData.str()<<") empty method");

    _needNewLHSJacobianLumped = false;
    // No implementation needed for interface.

    PYLITH_METHOD_END;
} // computeLHSJacobianLumpedInv


// ------------------------------------------------------------------------------------------------
// Compute weak form key part for face.
PetscInt
pylith::feassemble::_IntegratorInterface::computeWeakFormPart(const PetscInt part,
                                                              const PetscInt face) {
    return 3*part + face;
}


// ------------------------------------------------------------------------------------------------
// Compute residual.
void
pylith::feassemble::_IntegratorInterface::computeResidual(pylith::topology::Field* residual,
                                                          const pylith::feassemble::IntegratorInterface* integrator,
                                                          pylith::feassemble::Integrator::ResidualPart residualPart,
                                                          const pylith::problems::IntegrationData& integrationData) {
    PYLITH_METHOD_BEGIN;
    typedef InterfacePatches::keysmap_t keysmap_t;

    pythia::journal::debug_t debug(_IntegratorInterface::genericComponent);
    debug << pythia::journal::at(__HERE__)
          << "_IntegratorInterface::computeResidual(residual="<<typeid(residual).name()<<", integrator"<<typeid(integrator).name()
          <<", residualPart="<<residualPart<<", integrationData="<<integrationData.str()<<")"
          << pythia::journal::endl;

    assert(integrator);
    assert(residual);

    const pylith::topology::Field* solution = integrationData.getField(pylith::problems::IntegrationData::solution);
    assert(solution);
    const PylithReal t = integrationData.getScalar(pylith::problems::IntegrationData::time);
    const PylithReal dt = integrationData.getScalar(pylith::problems::IntegrationData::time_step);
    PetscVec solutionDotVec = NULL;
    if (residualPart == pylith::feassemble::Integrator::RESIDUAL_LHS) {
        const pylith::topology::Field* solutionDot = integrationData.getField(pylith::problems::IntegrationData::solution_dot);
        assert(solutionDot);
        solutionDotVec = solutionDot->getLocalVector();
    } // if

    integrator->_setKernelConstants(*solution, dt);

    // Loop over integration patches.
    PetscErrorCode err = 0;
    PetscDM dmSoln = solution->getDM();
    const InterfacePatches* patches = integrator->_integrationPatches;assert(patches);
    const keysmap_t& keysmap = patches->getKeys();
    for (keysmap_t::const_iterator iter = keysmap.begin(); iter != keysmap.end(); ++iter) {
        PetscFormKey weakFormKeys[3];
        weakFormKeys[0] = iter->second.negative.getPetscKey(*solution, residualPart);
        weakFormKeys[0].part = _IntegratorInterface::computeWeakFormPart(residualPart, IntegratorInterface::NEGATIVE_FACE);

        weakFormKeys[1] = iter->second.positive.getPetscKey(*solution, residualPart);
        weakFormKeys[1].part = _IntegratorInterface::computeWeakFormPart(residualPart, IntegratorInterface::POSITIVE_FACE);

        weakFormKeys[2] = iter->second.cohesive.getPetscKey(*solution, residualPart);
        weakFormKeys[2].part = _IntegratorInterface::computeWeakFormPart(residualPart, IntegratorInterface::COHESIVE_FACE);

        PetscIS patchCellsIS = NULL;
        PetscInt numPatchCells = 0;
        const PetscInt* patchCells = NULL;
        err = DMGetStratumIS(dmSoln, patches->getLabelName(), weakFormKeys[2].value, &patchCellsIS);PYLITH_CHECK_ERROR(err);
        err = ISGetSize(patchCellsIS, &numPatchCells);PYLITH_CHECK_ERROR(err);assert(numPatchCells > 0);
        err = ISGetIndices(patchCellsIS, &patchCells);PYLITH_CHECK_ERROR(err);assert(patchCells);
        assert(pylith::topology::MeshOps::isCohesiveCell(dmSoln, patchCells[0]));

        assert(solution->getLocalVector());
        assert(residual->getLocalVector());
        err = DMPlexComputeResidual_Hybrid_Internal(dmSoln, weakFormKeys, patchCellsIS, t, solution->getLocalVector(),
                                                    solutionDotVec, t, residual->getLocalVector(), NULL);PYLITH_CHECK_ERROR(err);
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
                                                          const pylith::problems::IntegrationData& integrationData) {
    PYLITH_METHOD_BEGIN;
    typedef InterfacePatches::keysmap_t keysmap_t;

    pythia::journal::debug_t debug(_IntegratorInterface::genericComponent);
    debug << pythia::journal::at(__HERE__)
          << "_IntegratorInterface::computeJacobian(jacobianMat="<<jacobianMat<<", precondMat"<<precondMat
          <<", integrator"<<typeid(integrator).name()<<", jacobianPart="<<jacobianPart<<", integrationData="<<integrationData.str()<<")"
          << pythia::journal::endl;

    assert(jacobianMat);
    assert(precondMat);
    assert(integrator);

    const pylith::topology::Field* solution = integrationData.getField(pylith::problems::IntegrationData::solution);
    assert(solution);
    const pylith::topology::Field* solutionDot = integrationData.getField(pylith::problems::IntegrationData::solution_dot);
    assert(solutionDot);
    const PylithReal t = integrationData.getScalar(pylith::problems::IntegrationData::time);
    const PylithReal dt = integrationData.getScalar(pylith::problems::IntegrationData::time_step);
    const PylithReal s_tshift = integrationData.getScalar(pylith::problems::IntegrationData::s_tshift);

    integrator->_setKernelConstants(*solution, dt);

    PetscErrorCode err;
    PetscDM dmSoln = solution->getDM();
    const InterfacePatches* patches = integrator->_integrationPatches;assert(patches);
    const keysmap_t& keysmap = patches->getKeys();
    for (keysmap_t::const_iterator iter = keysmap.begin(); iter != keysmap.end(); ++iter) {
        PetscFormKey weakFormKeys[3];
        weakFormKeys[0] = iter->second.negative.getPetscKey(*solution, jacobianPart);
        weakFormKeys[0].part = _IntegratorInterface::computeWeakFormPart(jacobianPart, IntegratorInterface::NEGATIVE_FACE);

        weakFormKeys[1] = iter->second.positive.getPetscKey(*solution, jacobianPart);
        weakFormKeys[1].part = _IntegratorInterface::computeWeakFormPart(jacobianPart, IntegratorInterface::POSITIVE_FACE);

        weakFormKeys[2] = iter->second.cohesive.getPetscKey(*solution, jacobianPart);
        weakFormKeys[2].part = _IntegratorInterface::computeWeakFormPart(jacobianPart, IntegratorInterface::COHESIVE_FACE);

        PetscIS patchCellsIS = NULL;
        PetscInt numPatchCells = 0;
        const PetscInt* patchCells = NULL;
        err = DMGetStratumIS(dmSoln, patches->getLabelName(), weakFormKeys[2].value, &patchCellsIS);PYLITH_CHECK_ERROR(err);
        err = ISGetSize(patchCellsIS, &numPatchCells);PYLITH_CHECK_ERROR(err);assert(numPatchCells > 0);
        err = ISGetIndices(patchCellsIS, &patchCells);PYLITH_CHECK_ERROR(err);assert(patchCells);
        assert(pylith::topology::MeshOps::isCohesiveCell(dmSoln, patchCells[0]));

        assert(solution->getLocalVector());
        err = DMPlexComputeJacobian_Hybrid_Internal(dmSoln, weakFormKeys, patchCellsIS, t, s_tshift, solution->getLocalVector(),
                                                    solutionDot->getLocalVector(), jacobianMat, precondMat,
                                                    NULL);PYLITH_CHECK_ERROR(err);
        err = ISRestoreIndices(patchCellsIS, &patchCells);PYLITH_CHECK_ERROR(err);
        err = ISDestroy(&patchCellsIS);PYLITH_CHECK_ERROR(err);
    }
    PYLITH_METHOD_END;
} // computeJacobian


// ------------------------------------------------------------------------------------------------
void
pylith::feassemble::_IntegratorInterface::addMaterialKernels(std::vector<ResidualKernels>* kernels,
                                                             const pylith::feassemble::IntegratorInterface::FaceEnum face,
                                                             const pylith::feassemble::FEKernelKey key,
                                                             const pylith::topology::Field& solution,
                                                             const std::vector<pylith::materials::Material*>& materials) {
    PYLITH_METHOD_BEGIN;

    const pylith::materials::Material* material =
        _IntegratorInterface::getMaterial(materials, key.getName(), key.getValue());assert(material);
    const std::vector<ResidualKernels>& kernelsMaterial = material->getInterfaceKernelsResidual(solution, face);
    kernels->insert(kernels->end(), kernelsMaterial.begin(), kernelsMaterial.end() );

    PYLITH_METHOD_END;
} // addMaterialKernels(Residual)


// ------------------------------------------------------------------------------------------------
void
pylith::feassemble::_IntegratorInterface::addMaterialKernels(std::vector<JacobianKernels>* kernels,
                                                             const pylith::feassemble::IntegratorInterface::FaceEnum face,
                                                             const pylith::feassemble::FEKernelKey key,
                                                             const pylith::topology::Field& solution,
                                                             const std::vector<pylith::materials::Material*>& materials) {
    PYLITH_METHOD_BEGIN;

    const pylith::materials::Material* material =
        _IntegratorInterface::getMaterial(materials, key.getName(), key.getValue());assert(material);
    const std::vector<JacobianKernels>& kernelsMaterial = material->getInterfaceKernelsJacobian(solution, face);
    kernels->insert(kernels->end(), kernelsMaterial.begin(), kernelsMaterial.end() );

    PYLITH_METHOD_END;
} // getMaterialKernels(Jacobian)


// ------------------------------------------------------------------------------------------------
pylith::materials::Material*
pylith::feassemble::_IntegratorInterface::getMaterial(const std::vector<pylith::materials::Material*> materials,
                                                      const char* labelName,
                                                      const int labelValue) {
    PYLITH_METHOD_BEGIN;

    const std::string& materialLabelName = pylith::topology::Mesh::getCellsLabelName();

    pylith::materials::Material* material = NULL;
    for (size_t i = 0; i < materials.size(); ++i) {
        if ((labelName == materialLabelName) && (labelValue == materials[i]->getMaterialId())) {
            material = materials[i];
            break;
        } // if
    } // for

    PYLITH_METHOD_RETURN(material);
} // getMaterial


// End of file
