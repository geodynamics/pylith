// =================================================================================================
// This code is part of PyLith, developed through the Computational Infrastructure
// for Geodynamics (https://github.com/geodynamics/pylith).
//
// Copyright (c) 2010-2025, University of California, Davis and the PyLith Development Team.
// All rights reserved.
//
// See https://mit-license.org/ and LICENSE.md and for license information.
// =================================================================================================

#include <portinfo>

#include "pylith/feassemble/IntegratorBoundary.hh" // implementation of object methods

#include "pylith/feassemble/DSLabelAccess.hh" // USES DSLabelAccess
#include "pylith/topology/Mesh.hh" // USES Mesh
#include "pylith/topology/MeshOps.hh" // USES createLowerDimMesh()
#include "pylith/topology/Field.hh" // USES Field
#include "pylith/problems/Physics.hh" // USES Physics
#include "pylith/feassemble/IntegrationData.hh" // USES IntegrationData
#include "pylith/topology/CoordsVisitor.hh" // USES CoordsVisitor::optimizeClosure()

#include "spatialdata/spatialdb/GravityField.hh" // HASA GravityField

#include "pylith/utils/journals.hh" // USES PYLITH_JOURNAL_*

#include "petscds.h" // USES PetscDS

#include <cassert> // USES assert()
#include <stdexcept> // USES std::runtime_error

// ------------------------------------------------------------------------------------------------
// Local "private" functions.
namespace pylith {
    namespace feassemble {
        class _IntegratorBoundary {
public:

            static const char* genericComponent;
        }; // _IntegratorBoundary
        const char* _IntegratorBoundary::genericComponent = "integratorboundary";

    } // feassemble
} // pylith

// ------------------------------------------------------------------------------------------------
// Default constructor.
pylith::feassemble::IntegratorBoundary::IntegratorBoundary(pylith::problems::Physics* const physics) :
    Integrator(physics),
    _boundaryMesh(NULL),
    _boundarySurfaceLabel(""),
    _subfieldName("") {
    GenericComponent::setName(_IntegratorBoundary::genericComponent);
} // constructor


// ------------------------------------------------------------------------------------------------
// Destructor.
pylith::feassemble::IntegratorBoundary::~IntegratorBoundary(void) {
    deallocate();
} // destructor


// ------------------------------------------------------------------------------------------------
// Deallocate PETSc and local data structures.
void
pylith::feassemble::IntegratorBoundary::deallocate(void) {
    PYLITH_METHOD_BEGIN;

    Integrator::deallocate();

    delete _boundaryMesh;_boundaryMesh = NULL;

    PYLITH_METHOD_END;
} // deallocate


// ------------------------------------------------------------------------------------------------
// Set name of solution subfield associated with boundary condition.
void
pylith::feassemble::IntegratorBoundary::setSubfieldName(const char* value) {
    PYLITH_JOURNAL_DEBUG("setSubfieldName(value="<<value<<")");

    if (!value || (0 == strlen(value))) {
        std::ostringstream msg;
        msg << "Empty string given for name of solution subfield for boundary condition '" << _boundarySurfaceLabel
            <<"'.";
        throw std::runtime_error(msg.str());
    } // if
    _subfieldName = value;
} // setSubfieldName


// ------------------------------------------------------------------------------------------------
// Get name of solution subfield associated with boundary condition.
const char*
pylith::feassemble::IntegratorBoundary::getSubfieldName(void) const {
    return _subfieldName.c_str();
} // getSubfieldName


// ------------------------------------------------------------------------------------------------
// Get mesh associated with integration domain.
const pylith::topology::Mesh&
pylith::feassemble::IntegratorBoundary::getPhysicsDomainMesh(void) const {
    assert(_boundaryMesh);
    return *_boundaryMesh;
} // domainMesh


// ------------------------------------------------------------------------------------------------
void
pylith::feassemble::IntegratorBoundary::setKernelsResidual(const std::vector<ResidualKernels>& kernels,
                                                           const pylith::topology::Field& solution) {
    PYLITH_METHOD_BEGIN;
    PYLITH_JOURNAL_DEBUG(_labelName<<"="<<_labelValue<<" setKernelsResidual(# kernels="<<kernels.size()<<")");

    PetscErrorCode err;
    DSLabelAccess dsLabel(solution.getDM(), _labelName.c_str(), _labelValue);
    for (size_t i = 0; i < kernels.size(); ++i) {
        const PetscInt i_field = solution.getSubfieldInfo(kernels[i].subfield.c_str()).index;
        const PetscInt i_part = kernels[i].part;
        if (dsLabel.weakForm()) {
            err = PetscWeakFormAddBdResidual(dsLabel.weakForm(), dsLabel.label(), dsLabel.value(), i_field, i_part,
                                             kernels[i].r0, kernels[i].r1);PYLITH_CHECK_ERROR(err);
        } // if

        switch (kernels[i].part) {
        case LHS:
            _hasLHSResidual = true;
            break;
        case RHS:
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


// ------------------------------------------------------------------------------------------------
void
pylith::feassemble::IntegratorBoundary::setKernelsDiagnosticField(const std::vector<ProjectKernels>& kernels) {
    PYLITH_METHOD_BEGIN;
    PYLITH_JOURNAL_DEBUG(_labelName<<"="<<_labelValue<<" setKernelsDiagnosticField(# kernels="<<kernels.size()<<")");

    _kernelsDiagnosticField = kernels;

    PYLITH_METHOD_END;
} // setKernelsDiagnosticField


// ------------------------------------------------------------------------------------------------
// Initialize integration domain, auxiliary field, and derived field. Update observers.
void
pylith::feassemble::IntegratorBoundary::initialize(const pylith::topology::Field& solution) {
    PYLITH_METHOD_BEGIN;
    PYLITH_JOURNAL_DEBUG(_labelName<<"="<<_labelValue<<" initialize(solution="<<solution.getLabel()<<")");

    delete _boundaryMesh;_boundaryMesh = pylith::topology::MeshOps::createLowerDimMesh(solution.getMesh(), _labelName.c_str(), _labelValue);
    assert(_boundaryMesh);
    pylith::topology::CoordsVisitor::optimizeClosure(_boundaryMesh->getDM());

    Integrator::initialize(solution);

    assert(_auxiliaryField);
    PetscErrorCode err;
    PetscDM dmSoln = solution.getDM();assert(dmSoln);
    PetscDMLabel dmLabel = NULL;
    err = DMGetLabel(dmSoln, _labelName.c_str(), &dmLabel);PYLITH_CHECK_ERROR(err);assert(dmLabel);
    err = DMSetAuxiliaryVec(dmSoln, dmLabel, _labelValue, LHS, _auxiliaryField->getLocalVector());PYLITH_CHECK_ERROR(err);
    err = DMSetAuxiliaryVec(dmSoln, dmLabel, _labelValue, RHS, _auxiliaryField->getLocalVector());PYLITH_CHECK_ERROR(err);

    pythia::journal::debug_t debug(GenericComponent::getName());
    if (debug.state()) {
        PYLITH_JOURNAL_DEBUG("Viewing auxiliary field.");
        _auxiliaryField->view("Auxiliary field");
    } // if

    PYLITH_METHOD_END;
} // initialize


// ---------------------------------------------------------------------------------------------------------------------
// Set auxiliary field values for current time.
void
pylith::feassemble::IntegratorBoundary::setState(const double t) {
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
pylith::feassemble::IntegratorBoundary::computeRHSResidual(pylith::topology::Field* residual,
                                                           const pylith::feassemble::IntegrationData& integrationData) {
    PYLITH_METHOD_BEGIN;
    PYLITH_JOURNAL_DEBUG(_labelName<<"="<<_labelValue<<" computeRHSResidual(residual="<<residual<<", integrationData="<<integrationData.str()<<")");
    if (!_hasRHSResidual) { PYLITH_METHOD_END;}
    assert(residual);

    const pylith::topology::Field* solution = integrationData.getField(pylith::feassemble::IntegrationData::solution);
    assert(solution);
    const PylithReal t = integrationData.getScalar(pylith::feassemble::IntegrationData::time);
    const PylithReal dt = integrationData.getScalar(pylith::feassemble::IntegrationData::time_step);

    DSLabelAccess dsLabel(solution->getDM(), _labelName.c_str(), _labelValue);
    _setKernelConstants(*solution, dt);

    PetscFormKey key;
    key.label = dsLabel.label();
    key.value = dsLabel.value();
    key.field = solution->getSubfieldInfo(_subfieldName.c_str()).index;
    key.part = pylith::feassemble::Integrator::RHS;

    PetscErrorCode err;
    assert(solution->getLocalVector());
    assert(residual->getLocalVector());
    PetscVec solutionDotVec = NULL;
    err = DMPlexComputeBdResidualSingle(dsLabel.dm(), t, dsLabel.weakForm(), key, solution->getLocalVector(), solutionDotVec,
                                        residual->getLocalVector());PYLITH_CHECK_ERROR(err);

    PYLITH_METHOD_END;
} // computeRHSResidual


// ------------------------------------------------------------------------------------------------
// Compute LHS residual for F(t,s,\dot{s}).
void
pylith::feassemble::IntegratorBoundary::computeLHSResidual(pylith::topology::Field* residual,
                                                           const pylith::feassemble::IntegrationData& integrationData) {
    PYLITH_METHOD_BEGIN;
    PYLITH_JOURNAL_DEBUG(_labelName<<"="<<_labelValue<<" computeLHSResidual(residual="<<residual<<", integrationData="<<integrationData.str()<<")");
    if (!_hasLHSResidual) { PYLITH_METHOD_END;}
    assert(residual);

    const pylith::topology::Field* solution = integrationData.getField(pylith::feassemble::IntegrationData::solution);
    assert(solution);
    const pylith::topology::Field* solutionDot = integrationData.getField(pylith::feassemble::IntegrationData::solution_dot);
    assert(solutionDot);
    const PylithReal t = integrationData.getScalar(pylith::feassemble::IntegrationData::time);
    const PylithReal dt = integrationData.getScalar(pylith::feassemble::IntegrationData::time_step);

    DSLabelAccess dsLabel(solution->getDM(), _labelName.c_str(), _labelValue);
    _setKernelConstants(*solution, dt);

    PetscFormKey key;
    key.label = dsLabel.label();
    key.value = dsLabel.value();
    key.field = solution->getSubfieldInfo(_subfieldName.c_str()).index;
    key.part = pylith::feassemble::Integrator::LHS;

    PetscErrorCode err;
    assert(solution->getLocalVector());
    assert(residual->getLocalVector());
    err = DMPlexComputeBdResidualSingle(dsLabel.dm(), t, dsLabel.weakForm(), key, solution->getLocalVector(),
                                        solutionDot->getLocalVector(), residual->getLocalVector());PYLITH_CHECK_ERROR(err);

    PYLITH_METHOD_END;
} // computeLHSResidual


// ------------------------------------------------------------------------------------------------
// Compute LHS Jacobian for F(t,s,\dot{s}).
void
pylith::feassemble::IntegratorBoundary::computeLHSJacobian(PetscMat jacobianMat,
                                                           PetscMat precondMat,
                                                           const pylith::feassemble::IntegrationData& integrationData) {
    PYLITH_METHOD_BEGIN;
    PYLITH_JOURNAL_DEBUG(_labelName<<"="<<_labelValue<<" computeLHSJacobian(jacobianMat="<<jacobianMat<<", precondMat="<<precondMat<<", integrationData="<<integrationData.str()<<") empty method");

    _needNewLHSJacobian = false;
    // No implementation needed for boundary.

    PYLITH_METHOD_END;
} // computeLHSJacobian


// ------------------------------------------------------------------------------------------------
// Compute LHS Jacobian for F(t,s,\dot{s}).
void
pylith::feassemble::IntegratorBoundary::computeLHSJacobianLumpedInv(pylith::topology::Field* jacobianInv,
                                                                    const pylith::feassemble::IntegrationData& integrationData) {
    PYLITH_METHOD_BEGIN;
    PYLITH_JOURNAL_DEBUG(_labelName<<"="<<_labelValue<<" computeLHSJacobianLumpedInv(jacobianInv="<<jacobianInv<<", integrationData="<<integrationData.str()<<") empty method");

    _needNewLHSJacobianLumped = false;
    // No implementation needed for boundary.

    PYLITH_METHOD_END;
} // computeLHSJacobianLumpedInv


// ------------------------------------------------------------------------------------------------
// Compute diagnostic field from auxiliary field.
void
pylith::feassemble::IntegratorBoundary::_computeDiagnosticField(void) {
    PYLITH_METHOD_BEGIN;
    PYLITH_JOURNAL_DEBUG("_computeDiagnosticField()");

    // We create the diagnostic field for all BCs, but we don't set the kernel for
    // BCs associated with analytical user functions.
    if (!_diagnosticField || !_auxiliaryField || !_kernelsDiagnosticField.size()) {
        PYLITH_METHOD_END;
    } // if

    assert(_auxiliaryField);
    assert(_diagnosticField);
    const PylithScalar t = 0.0;
    const PylithScalar dt = 0.0;
    _setKernelConstants(*_auxiliaryField, dt);

    const size_t numKernels = _kernelsDiagnosticField.size();
    assert(numKernels > 0);
    PetscBdPointFn** kernelsArray = (numKernels > 0) ? new PetscBdPointFn*[numKernels] : NULL;
    for (size_t iKernel = 0; iKernel < numKernels; ++iKernel) {
        const pylith::topology::Field::SubfieldInfo& sinfo = _diagnosticField->getSubfieldInfo(_kernelsDiagnosticField[iKernel].subfield.c_str());
        kernelsArray[sinfo.index] = _kernelsDiagnosticField[iKernel].f;
    } // for

    PetscErrorCode err = PETSC_SUCCESS;
    PetscDM diagnosticDM = _diagnosticField->getDM();
    PetscDMLabel diagnosticFieldLabel = NULL;
    const PetscInt labelValue = 1;
    err = DMGetLabel(diagnosticDM, "output", &diagnosticFieldLabel);PYLITH_CHECK_ERROR(err);
    err = DMProjectBdFieldLabelLocal(diagnosticDM, t, diagnosticFieldLabel, 1, &labelValue, PETSC_DETERMINE, NULL, _auxiliaryField->getLocalVector(), kernelsArray, INSERT_VALUES, _diagnosticField->getLocalVector());PYLITH_CHECK_ERROR(err);
    delete[] kernelsArray;kernelsArray = NULL;

    pythia::journal::debug_t debug(GenericComponent::getName());
    if (debug.state()) {
        PYLITH_JOURNAL_DEBUG("Viewing diagnostic field.");
        _diagnosticField->view("Diagnostic field");
    } // if

    PYLITH_METHOD_END;
} // _computeDiagnosticField


// End of file
