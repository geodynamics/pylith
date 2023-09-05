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
// Copyright (c) 2010-2022 University of California, Davis
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
#include "pylith/feassemble/IntegrationData.hh" // USES IntegrationData
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

            /** Compute residual using current kernels.
             *
             * @param[out] residual Field for residual.
             * @param[in] integrator Integrator for boundary.
             * @param[in] equationPart Equation part to compute.
             * @param[in] integrationData Data needed to integrate governing equations.
             */
            static
            void computeResidual(pylith::topology::Field* residual,
                                 const pylith::feassemble::IntegratorInterface* integrator,
                                 pylith::feassemble::Integrator::EquationPart equationPart,
                                 const pylith::feassemble::IntegrationData& integrationData);

            /** Compute Jacobian using current kernels.
             *
             * @param[out] jacobianMat PETSc Mat with Jacobian sparse matrix.
             * @param[out] precondMat PETSc Mat with Jacobian preconditioning sparse matrix.
             * @param[in] integrator Integrator for boundary.
             * @param[in] equationPart Equation part to compute.
             * @param[in] integrationData Data needed to integrate governing equations.
             */
            static
            void computeJacobian(PetscMat jacobianMat,
                                 PetscMat precondMat,
                                 const pylith::feassemble::IntegratorInterface* integrator,
                                 pylith::feassemble::Integrator::EquationPart equationPart,
                                 const pylith::feassemble::IntegrationData& integrationData);

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

            static const PetscInt max_face_enums; ///< Maximum number of fault faces (negative, positive, fault).
            static const PetscInt num_face_enums; ///< Number of fault faces (negative, positive, fault).
            static const PetscInt max_parts; ///< Maximum number of equation parts.
            static const PetscInt max_patches; ///< Maximum number of patches per fault.

            static const char* genericComponent;

        }; // _IntegratorInterface
        const char* _IntegratorInterface::genericComponent = "integratorinterface";
        const PetscInt _IntegratorInterface::max_face_enums = 10;
        const PetscInt _IntegratorInterface::num_face_enums = 3;
        const PetscInt _IntegratorInterface::max_parts = 10;
        const PetscInt _IntegratorInterface::max_patches = 10;

    } // feassemble
} // pylith

// ------------------------------------------------------------------------------------------------
// Default constructor.
pylith::feassemble::IntegratorInterface::IntegratorInterface(pylith::problems::Physics* const physics) :
    Integrator(physics),
    _interfaceMesh(NULL),
    _surfaceLabelName(""),
    _integrationPatches(NULL),
    _weightingDM(NULL),
    _weightingVec(NULL),
    _hasLHSResidualWeighted(false),
    _hasLHSJacobianWeighted(false) {
    GenericComponent::setName(_IntegratorInterface::genericComponent);
    _labelValue = 100;
    _labelName = pylith::topology::Mesh::cells_label_name;
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
    DMDestroy(&_weightingDM);
    VecDestroy(&_weightingVec);

    PYLITH_METHOD_END;
} // deallocate


// ------------------------------------------------------------------------------------------------
// Set label marking boundary associated with boundary condition surface.
void
pylith::feassemble::IntegratorInterface::setSurfaceLabelName(const char* value) {
    PYLITH_JOURNAL_DEBUG("setSurfaceLabelName(value="<<value<<")");

    if (strlen(value) == 0) {
        throw std::runtime_error("Empty string given for boundary condition label.");
    } // if

    _surfaceLabelName = value;
} // setSurfaceLabelName


// ------------------------------------------------------------------------------------------------
// Get label marking boundary associated with boundary condition surface.
const char*
pylith::feassemble::IntegratorInterface::getSurfaceLabelName(void) const {
    return _surfaceLabelName.c_str();
} // getSurfaceLabelName


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

    if (patches) {
        const size_t numPatches = patches->getKeys().size();
        const size_t maxNumPatches = _IntegratorInterface::max_patches;
        if (numPatches > maxNumPatches) {
            std::ostringstream msg;
            msg << "Number of integration patches ("<< numPatches << ") for interface integration "
                << _labelName << "=" << _labelValue
                << " exceeds maximum number of allowed patches (" << maxNumPatches << ").\n"
                << "Consolidate bulk materials or adjust maximum number of patches in "
                << "pylith::feassemble::IntegratorInterface.";
            throw std::range_error(msg.str());
        } // if
    } // if

    delete _integrationPatches;_integrationPatches = patches;

    PYLITH_METHOD_END;
} // setIntegrationPatches


// ------------------------------------------------------------------------------------------------
// Get integration patches.
const pylith::feassemble::InterfacePatches*
pylith::feassemble::IntegratorInterface::getIntegrationPatches(void) const {
    return _integrationPatches;
} // getIntegrationPatches


// ------------------------------------------------------------------------------------------------
// Set kernels for residual.
void
pylith::feassemble::IntegratorInterface::setKernels(const std::vector<ResidualKernels>& kernels,
                                                    const pylith::topology::Field& solution,
                                                    const std::vector<pylith::materials::Material*>& materials) {
    PYLITH_METHOD_BEGIN;
    PYLITH_JOURNAL_DEBUG(_labelName<<"="<<_labelValue<<" setKernels(# kernels="<<kernels.size()<<")");
    typedef InterfacePatches::keysmap_t keysmap_t;

    PetscErrorCode err;
    const keysmap_t& keysmap = _integrationPatches->getKeys();
    for (keysmap_t::const_iterator iter = keysmap.begin(); iter != keysmap.end(); ++iter) {
        const PetscWeakForm weakForm = iter->second.cohesive.getWeakForm();
        const PetscInt patchValue = iter->second.cohesive.getValue();
        std::vector<ResidualKernels> patchKernels(kernels);

        _IntegratorInterface::addMaterialKernels(&patchKernels, pylith::feassemble::IntegratorInterface::NEGATIVE_FACE,
                                                 iter->second.negative, solution, materials);
        _IntegratorInterface::addMaterialKernels(&patchKernels, pylith::feassemble::IntegratorInterface::POSITIVE_FACE,
                                                 iter->second.positive, solution, materials);

        for (size_t i = 0; i < patchKernels.size(); ++i) {
            PetscFormKey key;

            const PetscInt interfacePart = getWeakFormPart(patchKernels[i].part, patchKernels[i].face, patchValue);

            switch (kernels[i].face) {
            case IntegratorInterface::NEGATIVE_FACE:
                key = iter->second.negative.getPetscKey(solution, interfacePart, patchKernels[i].subfield.c_str());
                break;
            case IntegratorInterface::POSITIVE_FACE:
                key = iter->second.positive.getPetscKey(solution, interfacePart, patchKernels[i].subfield.c_str());
                break;
            case IntegratorInterface::FAULT_FACE:
                key = iter->second.cohesive.getPetscKey(solution, interfacePart, patchKernels[i].subfield.c_str());
                break;
            default:
                PYLITH_JOURNAL_LOGICERROR("Unknown integration face ("<<patchKernels[i].face<<").");
            } // switch
            if (weakForm) {
                err = PetscWeakFormAddBdResidual(weakForm, key.label, key.value, key.field, key.part,
                                                 patchKernels[i].r0, patchKernels[i].r1);PYLITH_CHECK_ERROR(err);
            } // if

            switch (patchKernels[i].part) {
            case LHS:
                _hasLHSResidual = true;
                break;
            case LHS_WEIGHTED:
                _hasLHSResidualWeighted = true;
                break;
            case RHS:
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
} // setKernels


// ------------------------------------------------------------------------------------------------
// Set kernels for Jacobian.
void
pylith::feassemble::IntegratorInterface::setKernels(const std::vector<JacobianKernels>& kernels,
                                                    const pylith::topology::Field& solution,
                                                    const std::vector<pylith::materials::Material*>& materials) {
    PYLITH_METHOD_BEGIN;
    PYLITH_JOURNAL_DEBUG(_labelName<<"="<<_labelValue<<" setKernels(# kernels="<<kernels.size()<<")");
    typedef InterfacePatches::keysmap_t keysmap_t;

    PetscErrorCode err = 0;
    const keysmap_t& keysmap = _integrationPatches->getKeys();
    for (keysmap_t::const_iterator iter = keysmap.begin(); iter != keysmap.end(); ++iter) {
        const PetscWeakForm weakForm = iter->second.cohesive.getWeakForm();
        const PetscInt patchValue = iter->second.cohesive.getValue();
        std::vector<JacobianKernels> patchKernels(kernels);

        _IntegratorInterface::addMaterialKernels(&patchKernels, pylith::feassemble::IntegratorInterface::NEGATIVE_FACE,
                                                 iter->second.negative, solution, materials);
        _IntegratorInterface::addMaterialKernels(&patchKernels, pylith::feassemble::IntegratorInterface::POSITIVE_FACE,
                                                 iter->second.positive, solution, materials);

        PetscInt numFields = 0;
        err = PetscWeakFormGetNumFields(weakForm, &numFields);PYLITH_CHECK_ERROR(err);

        for (size_t i = 0; i < patchKernels.size(); ++i) {
            PetscFormKey key;
            const PetscInt interfacePart = getWeakFormPart(patchKernels[i].part, patchKernels[i].face, patchValue);

            switch (patchKernels[i].face) {
            case IntegratorInterface::NEGATIVE_FACE:
                key = iter->second.negative.getPetscKey(solution, interfacePart, patchKernels[i].subfieldTrial.c_str(),
                                                        patchKernels[i].subfieldBasis.c_str());
                break;
            case IntegratorInterface::POSITIVE_FACE:
                key = iter->second.positive.getPetscKey(solution, interfacePart, patchKernels[i].subfieldTrial.c_str(),
                                                        patchKernels[i].subfieldBasis.c_str());
                break;
            case IntegratorInterface::FAULT_FACE:
                key = iter->second.cohesive.getPetscKey(solution, interfacePart, patchKernels[i].subfieldTrial.c_str(),
                                                        patchKernels[i].subfieldBasis.c_str());
                break;
            default:
                PYLITH_JOURNAL_LOGICERROR("Unknown integration face.");
            } // switch
            const PetscInt i_trial = key.field / numFields;
            const PetscInt i_basis = key.field % numFields;
            if (weakForm) {
                err = PetscWeakFormAddBdJacobian(weakForm, key.label, key.value, i_trial, i_basis, key.part,
                                                 patchKernels[i].j0, patchKernels[i].j1, patchKernels[i].j2, patchKernels[i].j3);PYLITH_CHECK_ERROR(err);
            } // if

            switch (patchKernels[i].part) {
            case LHS:
                _hasLHSJacobian = true;
                break;
            case LHS_WEIGHTED:
                _hasLHSJacobianWeighted = true;
                break;
            case LHS_LUMPED_INV:
                _hasLHSJacobianLumped = true;
                break;
            default:
                PYLITH_JOURNAL_LOGICERROR("Unknown Jacobian part " << patchKernels[i].part <<".");
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
void
pylith::feassemble::IntegratorInterface::setKernelsUpdateStateVars(const std::vector<ProjectKernels>& kernels) {
    PYLITH_METHOD_BEGIN;
    PYLITH_JOURNAL_DEBUG(_labelName<<"="<<_labelValue<<" setKernelsUpdateStateVars(# kernels="<<kernels.size()<<")");

    _kernelsUpdateStateVars = kernels;

    PYLITH_METHOD_END;
} // setKernelsUpdateStateVars


// ------------------------------------------------------------------------------------------------
void
pylith::feassemble::IntegratorInterface::setKernelsDerivedField(const std::vector<ProjectKernels>& kernels) {
    PYLITH_METHOD_BEGIN;
    PYLITH_JOURNAL_DEBUG(_labelName<<"="<<_labelValue<<" setKernelsDerivedField(# kernels="<<kernels.size()<<")");

    _kernelsDerivedField = kernels;

    PYLITH_METHOD_END;
} // setKernelsDerivedField


// ------------------------------------------------------------------------------------------------
// Compute weak form key part for face.
PetscInt
pylith::feassemble::IntegratorInterface::getWeakFormPart(const PetscInt part,
                                                         const PetscInt face,
                                                         const PetscInt patch) const {
    const PetscInt max_face_enums = _IntegratorInterface::max_face_enums;
    const PetscInt num_face_enums = _IntegratorInterface::num_face_enums;
    const PetscInt max_parts = _IntegratorInterface::max_parts;
    const PetscInt max_patches = _IntegratorInterface::max_patches;
    return _labelValue*(max_parts*max_face_enums*max_patches) + part*num_face_enums*max_patches + face*max_patches + patch;
} // getWeakFormPart


// ------------------------------------------------------------------------------------------------
// Initialize integration domain, auxiliary field, and derived field. Update observers.
void
pylith::feassemble::IntegratorInterface::initialize(const pylith::topology::Field& solution) {
    PYLITH_METHOD_BEGIN;
    PYLITH_JOURNAL_DEBUG(_labelName<<"="<<_labelValue<<" intialize(solution="<<solution.getLabel()<<")");

    const bool isSubmesh = true;
    delete _interfaceMesh;_interfaceMesh = new pylith::topology::Mesh(isSubmesh);assert(_interfaceMesh);
    pylith::faults::TopologyOps::createFaultParallel(_interfaceMesh, solution.getMesh(), _labelValue, _labelName.c_str(),
                                                     _surfaceLabelName.c_str());
    pylith::topology::MeshOps::checkTopology(*_interfaceMesh);
    pylith::topology::CoordsVisitor::optimizeClosure(_interfaceMesh->getDM());

    Integrator::initialize(solution);

    PetscErrorCode err = 0;
    PetscDM dmSoln = solution.getDM();
    typedef InterfacePatches::keysmap_t keysmap_t;
    const keysmap_t& keysmap = _integrationPatches->getKeys();
    const pylith::topology::Field* auxiliaryField = getAuxiliaryField();assert(auxiliaryField);
    for (keysmap_t::const_iterator iter = keysmap.begin(); iter != keysmap.end(); ++iter) {
        const PetscInt patchValue = iter->second.cohesive.getValue();
        const size_t numParts = 3;
        const EquationPart equationParts[numParts] = {
            pylith::feassemble::Integrator::RHS,
            pylith::feassemble::Integrator::LHS,
            pylith::feassemble::Integrator::LHS_WEIGHTED,
        };
        for (size_t i = 0; i < numParts; ++i) {
            PetscFormKey key = iter->second.cohesive.getPetscKey(solution, equationParts[i]);
            PetscInt part = getWeakFormPart(key.part, IntegratorInterface::FAULT_FACE, patchValue);
            err = DMSetAuxiliaryVec(dmSoln, key.label, key.value, part,
                                    auxiliaryField->getLocalVector());PYLITH_CHECK_ERROR(err);
        } // for
    } // for

    PYLITH_METHOD_END;
} // initialize


// ------------------------------------------------------------------------------------------------
// Set auxiliary field values for current time.
void
pylith::feassemble::IntegratorInterface::setState(const double t) {
    PYLITH_METHOD_BEGIN;
    PYLITH_JOURNAL_DEBUG(_labelName<<"="<<_labelValue<<" setState(t="<<t<<")");

    Integrator::setState(t);

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
} // setState


// ------------------------------------------------------------------------------------------------
// Compute RHS residual for G(t,s).
void
pylith::feassemble::IntegratorInterface::computeRHSResidual(pylith::topology::Field* residual,
                                                            const pylith::feassemble::IntegrationData& integrationData) {
    PYLITH_METHOD_BEGIN;
    PYLITH_JOURNAL_DEBUG(_labelName<<"="<<_labelValue<<" computeRHSResidual(residual="<<residual<<", integrationData="<<integrationData.str()<<")");
    if (!_hasRHSResidual) { PYLITH_METHOD_END;}

    _IntegratorInterface::computeResidual(residual, this, pylith::feassemble::Integrator::RHS, integrationData);

    PYLITH_METHOD_END;
} // computeRHSResidual


// ------------------------------------------------------------------------------------------------
// Compute LHS residual for F(t,s,\dot{s}).
void
pylith::feassemble::IntegratorInterface::computeLHSResidual(pylith::topology::Field* residual,
                                                            const pylith::feassemble::IntegrationData& integrationData) {
    PYLITH_METHOD_BEGIN;
    PYLITH_JOURNAL_DEBUG(_labelName<<"="<<_labelValue<<" computeLHSResidual(residual="<<residual<<", integrationData="<<integrationData.str()<<")");
    if (!_hasLHSResidual) { PYLITH_METHOD_END;}

    if (_hasLHSResidual) {
        const pylith::feassemble::Integrator::EquationPart equationPart = pylith::feassemble::Integrator::LHS;
        _IntegratorInterface::computeResidual(residual, this, equationPart, integrationData);
    } // if

    if (_hasLHSResidualWeighted) {
        { // KLUDGE
            PetscErrorCode err = 0;
            const pylith::topology::Field* daeWeighting =
                integrationData.getField(pylith::feassemble::IntegrationData::dae_mass_weighting);
            const pylith::topology::Field* solution =
                integrationData.getField(pylith::feassemble::IntegrationData::solution);
            PetscDM dmSoln = solution->getDM();
            typedef InterfacePatches::keysmap_t keysmap_t;
            const keysmap_t& keysmap = _integrationPatches->getKeys();

            for (keysmap_t::const_iterator iter = keysmap.begin(); iter != keysmap.end(); ++iter) {
                const PetscInt patchValue = iter->second.cohesive.getValue();
                const PetscFormKey key = iter->second.cohesive.getPetscKey(*solution, LHS_WEIGHTED);
                PetscInt part = getWeakFormPart(key.part, IntegratorInterface::FAULT_FACE, patchValue);
                err = DMSetAuxiliaryVec(dmSoln, key.label, -key.value, part, daeWeighting->getLocalVector());PYLITH_CHECK_ERROR(err);

                part = getWeakFormPart(key.part, IntegratorInterface::NEGATIVE_FACE, patchValue);
                err = DMSetAuxiliaryVec(dmSoln, key.label, -key.value, part, daeWeighting->getLocalVector());PYLITH_CHECK_ERROR(err);

                part = getWeakFormPart(key.part, IntegratorInterface::POSITIVE_FACE, patchValue);
                err = DMSetAuxiliaryVec(dmSoln, key.label, -key.value, part, daeWeighting->getLocalVector());PYLITH_CHECK_ERROR(err);
            } // for
        } // KLUDGE

        const pylith::feassemble::Integrator::EquationPart equationPart = pylith::feassemble::Integrator::LHS_WEIGHTED;
        _IntegratorInterface::computeResidual(residual, this, equationPart, integrationData);
    } // if

    PYLITH_METHOD_END;
} // computeLHSResidual


// ------------------------------------------------------------------------------------------------
// Compute LHS Jacobian for F(t,s,\dot{s}).
void
pylith::feassemble::IntegratorInterface::computeLHSJacobian(PetscMat jacobianMat,
                                                            PetscMat precondMat,
                                                            const pylith::feassemble::IntegrationData& integrationData) {
    PYLITH_METHOD_BEGIN;
    PYLITH_JOURNAL_DEBUG(_labelName<<"="<<_labelValue<<" computeLHSJacobian(jacobianMat="<<jacobianMat<<", precondMat="<<precondMat<<", integrationData="<<integrationData.str()<<")");

    _needNewLHSJacobian = false;

    if (_hasLHSJacobian) {
        pylith::feassemble::Integrator::EquationPart equationPart = pylith::feassemble::Integrator::LHS;
        _IntegratorInterface::computeJacobian(jacobianMat, precondMat, this, equationPart, integrationData);
    } // if

    if (_hasLHSJacobianWeighted) {
        pylith::feassemble::Integrator::EquationPart equationPart = pylith::feassemble::Integrator::LHS_WEIGHTED;
        _IntegratorInterface::computeJacobian(jacobianMat, precondMat, this, equationPart, integrationData);
    } // if

    PYLITH_METHOD_END;
} // computeLHSJacobian


// ------------------------------------------------------------------------------------------------
// Compute LHS Jacobian for F(t,s,\dot{s}).
void
pylith::feassemble::IntegratorInterface::computeLHSJacobianLumpedInv(pylith::topology::Field* jacobianInv,
                                                                     const pylith::feassemble::IntegrationData& integrationData) {
    PYLITH_METHOD_BEGIN;
    PYLITH_JOURNAL_DEBUG(_labelName<<"="<<_labelValue<<" computeLHSJacobianLumpedInv(jacobianInv="<<jacobianInv<<", integrationData="<<integrationData.str()<<") empty method");

    _needNewLHSJacobianLumped = false;
    // No implementation needed for interface.

    PYLITH_METHOD_END;
} // computeLHSJacobianLumpedInv


// ------------------------------------------------------------------------------------------------
// Compute residual.
void
pylith::feassemble::_IntegratorInterface::computeResidual(pylith::topology::Field* residual,
                                                          const pylith::feassemble::IntegratorInterface* integrator,
                                                          pylith::feassemble::Integrator::EquationPart equationPart,
                                                          const pylith::feassemble::IntegrationData& integrationData) {
    PYLITH_METHOD_BEGIN;
    typedef InterfacePatches::keysmap_t keysmap_t;

    pythia::journal::debug_t debug(_IntegratorInterface::genericComponent);
    debug << pythia::journal::at(__HERE__)
          << "_IntegratorInterface::computeResidual(residual="<<typeid(residual).name()<<", integrator"<<typeid(integrator).name()
          <<", equationPart="<<equationPart<<", integrationData="<<integrationData.str()<<")"
          << pythia::journal::endl;

    assert(integrator);
    assert(residual);

    const pylith::topology::Field* solution = integrationData.getField(pylith::feassemble::IntegrationData::solution);
    assert(solution);
    const PylithReal t = integrationData.getScalar(pylith::feassemble::IntegrationData::time);
    const PylithReal dt = integrationData.getScalar(pylith::feassemble::IntegrationData::time_step);
    PetscVec solutionDotVec = NULL;
    if ((equationPart == pylith::feassemble::Integrator::LHS) ||
        (equationPart == pylith::feassemble::Integrator::LHS_WEIGHTED) ) {
        const pylith::topology::Field* solutionDot = integrationData.getField(pylith::feassemble::IntegrationData::solution_dot);
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
        const PetscInt patchValue = iter->second.cohesive.getValue();

        PetscFormKey weakFormKeys[3];
        weakFormKeys[0] = iter->second.negative.getPetscKey(*solution, equationPart);
        weakFormKeys[0].part = integrator->getWeakFormPart(equationPart, IntegratorInterface::NEGATIVE_FACE, patchValue);

        weakFormKeys[1] = iter->second.positive.getPetscKey(*solution, equationPart);
        weakFormKeys[1].part = integrator->getWeakFormPart(equationPart, IntegratorInterface::POSITIVE_FACE, patchValue);

        weakFormKeys[2] = iter->second.cohesive.getPetscKey(*solution, equationPart);
        weakFormKeys[2].part = integrator->getWeakFormPart(equationPart, IntegratorInterface::FAULT_FACE, patchValue);

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
                                                          pylith::feassemble::Integrator::EquationPart equationPart,
                                                          const pylith::feassemble::IntegrationData& integrationData) {
    PYLITH_METHOD_BEGIN;
    typedef InterfacePatches::keysmap_t keysmap_t;

    pythia::journal::debug_t debug(_IntegratorInterface::genericComponent);
    debug << pythia::journal::at(__HERE__)
          << "_IntegratorInterface::computeJacobian(jacobianMat="<<jacobianMat<<", precondMat"<<precondMat
          <<", integrator"<<typeid(integrator).name()<<", equationPart="<<equationPart<<", integrationData="<<integrationData.str()<<")"
          << pythia::journal::endl;

    assert(jacobianMat);
    assert(precondMat);
    assert(integrator);

    const pylith::topology::Field* solution = integrationData.getField(pylith::feassemble::IntegrationData::solution);
    assert(solution);
    const pylith::topology::Field* solutionDot = integrationData.getField(pylith::feassemble::IntegrationData::solution_dot);
    assert(solutionDot);
    const PylithReal t = integrationData.getScalar(pylith::feassemble::IntegrationData::time);
    const PylithReal dt = integrationData.getScalar(pylith::feassemble::IntegrationData::time_step);
    const PylithReal s_tshift = integrationData.getScalar(pylith::feassemble::IntegrationData::s_tshift);

    integrator->_setKernelConstants(*solution, dt);

    PetscErrorCode err;
    PetscDM dmSoln = solution->getDM();
    const InterfacePatches* patches = integrator->_integrationPatches;assert(patches);
    const keysmap_t& keysmap = patches->getKeys();
    for (keysmap_t::const_iterator iter = keysmap.begin(); iter != keysmap.end(); ++iter) {
        const PetscInt patchValue = iter->second.cohesive.getValue();

        PetscFormKey weakFormKeys[3];
        weakFormKeys[0] = iter->second.negative.getPetscKey(*solution, equationPart);
        weakFormKeys[0].part = integrator->getWeakFormPart(equationPart, IntegratorInterface::NEGATIVE_FACE, patchValue);

        weakFormKeys[1] = iter->second.positive.getPetscKey(*solution, equationPart);
        weakFormKeys[1].part = integrator->getWeakFormPart(equationPart, IntegratorInterface::POSITIVE_FACE, patchValue);

        weakFormKeys[2] = iter->second.cohesive.getPetscKey(*solution, equationPart);
        weakFormKeys[2].part = integrator->getWeakFormPart(equationPart, IntegratorInterface::FAULT_FACE, patchValue);

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
// Update state variables as needed.
void
pylith::feassemble::IntegratorInterface::_updateStateVars(const PylithReal t,
                                                          const PylithReal dt,
                                                          const pylith::topology::Field& solution) {
    PYLITH_METHOD_BEGIN;
    PYLITH_JOURNAL_DEBUG("_updateStateVars(t="<<t<<", dt="<<dt<<", solution="<<solution.getLabel()<<")");

    if (0 == _kernelsUpdateStateVars.size()) {
        PYLITH_METHOD_END;
    } // if

    PYLITH_JOURNAL_LOGICERROR("_updateStateVars() not implemented.");

    PYLITH_METHOD_END;
} // _updateStateVars


// ------------------------------------------------------------------------------------------------
// Compute field derived from solution and auxiliary field.
void
pylith::feassemble::IntegratorInterface::_computeDerivedField(const PylithReal t,
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
    PetscBdPointFunc* kernelsArray = (numKernels > 0) ? new PetscBdPointFunc[numKernels] : NULL;
    for (size_t iKernel = 0; iKernel < numKernels; ++iKernel) {
        const pylith::topology::Field::SubfieldInfo& sinfo = _derivedField->getSubfieldInfo(_kernelsDerivedField[iKernel].subfield.c_str());
        kernelsArray[sinfo.index] = _kernelsDerivedField[iKernel].f;
    } // for

    PetscErrorCode err = 0;

    PetscDM derivedDM = _derivedField->getDM();
    assert(_auxiliaryField);
    PetscDMLabel dmLabel = NULL;
    PetscInt labelValue = 0;
    const PetscInt part = 0;
    err = DMSetAuxiliaryVec(derivedDM, dmLabel, labelValue, part, _auxiliaryField->getLocalVector());PYLITH_CHECK_ERROR(err);

    PetscDMLabel derivedFieldLabel = NULL;
    err = DMGetLabel(derivedDM, "output", &derivedFieldLabel);PYLITH_CHECK_ERROR(err);
    labelValue = 1;
    err = DMProjectBdFieldLabelLocal(derivedDM, t, derivedFieldLabel, 1, &labelValue, PETSC_DETERMINE, NULL, solution.getLocalVector(), kernelsArray, INSERT_VALUES, _derivedField->getLocalVector());PYLITH_CHECK_ERROR(err);
    delete[] kernelsArray;kernelsArray = NULL;

    pythia::journal::debug_t debug(GenericComponent::getName());
    if (debug.state()) {
        PYLITH_JOURNAL_DEBUG("Viewing derived field.");
        _derivedField->view("Derived field");
    } // if

    PYLITH_METHOD_END;
} // _computeDerivedField


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

    pylith::materials::Material* material = NULL;
    for (size_t i = 0; i < materials.size(); ++i) {
        assert(materials[i]);
        if ((std::string(labelName) == std::string(materials[i]->getLabelName())) && (labelValue == materials[i]->getLabelValue())) {
            material = materials[i];
            break;
        } // if
    } // for

    PYLITH_METHOD_RETURN(material);
} // getMaterial


// End of file
