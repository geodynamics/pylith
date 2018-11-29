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

#include <stdexcept> // USES std::runtime_error

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
             * @param[in] labelValue Value of label for interface cells for kernels;
             * @param[in] t Current time.
             * @param[in] dt Current time step.
             * @param[in] solution Field with current trial solution.
             * @param[in] solutionDot Field with time derivative of current trial solution.
             */
            static
            void computeResidual(pylith::topology::Field* residual,
                                 const pylith::feassemble::IntegratorInterface* integrator,
                                 const std::vector<pylith::feassemble::IntegratorInterface::ResidualKernels>& kernels,
                                 const PylithInt labelValue,
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
             * @param[in] labelValue Value of label for interface cells for kernels;
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
                                 const PylithInt labelValue,
                                 const PylithReal t,
                                 const PylithReal dt,
                                 const PylithReal s_tshift,
                                 const pylith::topology::Field& solution,
                                 const pylith::topology::Field& solutionDot);

        }; // _IntegratorInterface

    } // feassemble
} // pylith

// ---------------------------------------------------------------------------------------------------------------------
// Default constructor.
pylith::feassemble::IntegratorInterface::IntegratorInterface(pylith::problems::Physics* const physics) :
    Integrator(physics),
    _interfaceMesh(NULL),
    _interfaceId(100),
    _interfaceLabel("") {}


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
// Set value of label material-id used to identify interface cells.
void
pylith::feassemble::IntegratorInterface::setInterfaceId(const int value) {
    _interfaceId = value;
} // setInterfaceId


// ---------------------------------------------------------------------------------------------------------------------
// Get value of label material-id used to identify interface cells.
int
pylith::feassemble::IntegratorInterface::getInterfaceId(void) const {
    return _interfaceId;
} // getInterfaceId


// ---------------------------------------------------------------------------------------------------------------------
// Set label marking boundary associated with boundary condition surface.
void
pylith::feassemble::IntegratorInterface::setSurfaceMarkerLabel(const char* value) {
    if (strlen(value) == 0) {
        throw std::runtime_error("Empty string given for boundary condition label.");
    } // if

    _interfaceLabel = value;
} // setSurfaceMarkerLabel


// ---------------------------------------------------------------------------------------------------------------------
// Get label marking boundary associated with boundary condition surface.
const char*
pylith::feassemble::IntegratorInterface::getSurfaceMarkerLabel(void) const {
    return _interfaceLabel.c_str();
} // getSurfaceMarkerLabel


// ---------------------------------------------------------------------------------------------------------------------
// Get mesh associated with integration domain.
const pylith::topology::Mesh&
pylith::feassemble::IntegratorInterface::getIntegrationDomainMesh(void) const {
    assert(_interfaceMesh);
    return *_interfaceMesh;
} // domainMesh


// ---------------------------------------------------------------------------------------------------------------------
// Set kernels for RHS residual for the positive side of the interface.
void
pylith::feassemble::IntegratorInterface::setKernelsRHSResidualPos(const std::vector<ResidualKernels>& kernels) {
    PYLITH_METHOD_BEGIN;
    PYLITH_JOURNAL_DEBUG("setKernelsRHSResidual(# kernels="<<kernels.size()<<")");

    _kernelsRHSResidualPos = kernels;

    PYLITH_METHOD_END;
} // setKernelsRHSResidual


// ---------------------------------------------------------------------------------------------------------------------
// Set kernels for RHS residual for the negative side of the interface.
void
pylith::feassemble::IntegratorInterface::setKernelsRHSResidualNeg(const std::vector<ResidualKernels>& kernels) {
    PYLITH_METHOD_BEGIN;
    PYLITH_JOURNAL_DEBUG("setKernelsRHSResidualNeg(# kernels="<<kernels.size()<<")");

    _kernelsRHSResidualNeg = kernels;

    PYLITH_METHOD_END;
} // setKernelsRHSResidualNeg


// ---------------------------------------------------------------------------------------------------------------------
// Set kernels for RHS Jacobian for the positive side of the interface.
void
pylith::feassemble::IntegratorInterface::setKernelsRHSJacobianPos(const std::vector<JacobianKernels>& kernels) {
    PYLITH_METHOD_BEGIN;
    PYLITH_JOURNAL_DEBUG("setKernelsRHSJacobianPos(# kernels="<<kernels.size()<<")");

    _kernelsRHSJacobianPos = kernels;

    PYLITH_METHOD_END;
} // setKernelsRHSJacobianPOs


// ---------------------------------------------------------------------------------------------------------------------
// Set kernels for RHS Jacobian for the negative side of the interface.
void
pylith::feassemble::IntegratorInterface::setKernelsRHSJacobianNeg(const std::vector<JacobianKernels>& kernels) {
    PYLITH_METHOD_BEGIN;
    PYLITH_JOURNAL_DEBUG("setKernelsRHSJacobianNeg(# kernels="<<kernels.size()<<")");

    _kernelsRHSJacobianNeg = kernels;

    PYLITH_METHOD_END;
} // setKernelsRHSJacobianNegative


// ---------------------------------------------------------------------------------------------------------------------
// Set kernels for LHS residual for the positive side of the interface.
void
pylith::feassemble::IntegratorInterface::setKernelsLHSResidualPos(const std::vector<ResidualKernels>& kernels) {
    PYLITH_METHOD_BEGIN;
    PYLITH_JOURNAL_DEBUG("setKernelsLHSResidualPos(# kernels="<<kernels.size()<<")");

    _kernelsLHSResidualPos = kernels;

    PYLITH_METHOD_END;
} // setKernelsLHSResidualPos


// ---------------------------------------------------------------------------------------------------------------------
// Set kernels for LHS residual for the negative side of the interface.
void
pylith::feassemble::IntegratorInterface::setKernelsLHSResidualNeg(const std::vector<ResidualKernels>& kernels) {
    PYLITH_METHOD_BEGIN;
    PYLITH_JOURNAL_DEBUG("setKernelsLHSResidualNeg(# kernels="<<kernels.size()<<")");

    _kernelsLHSResidualNeg = kernels;

    PYLITH_METHOD_END;
} // setKernelsLHSResidualNeg


// ---------------------------------------------------------------------------------------------------------------------
// Set kernels for LHS Jacobian for the positive side of the interface.
void
pylith::feassemble::IntegratorInterface::setKernelsLHSJacobianPos(const std::vector<JacobianKernels>& kernels) {
    PYLITH_METHOD_BEGIN;
    PYLITH_JOURNAL_DEBUG("setKernelsLHSJacobianPos(# kernels="<<kernels.size()<<")");

    _kernelsLHSJacobianPos = kernels;

    PYLITH_METHOD_END;
} // setKernelsLHSJacobianPos


// ---------------------------------------------------------------------------------------------------------------------
// Set kernels for LHS Jacobian for the negative side of the interface.
void
pylith::feassemble::IntegratorInterface::setKernelsLHSJacobianNeg(const std::vector<JacobianKernels>& kernels) {
    PYLITH_METHOD_BEGIN;
    PYLITH_JOURNAL_DEBUG("setKernelsLHSJacobianNeg(# kernels="<<kernels.size()<<")");

    _kernelsLHSJacobianNeg = kernels;

    PYLITH_METHOD_END;
} // setKernelsLHSJacobianNeg


// ---------------------------------------------------------------------------------------------------------------------
// Initialize integration domain, auxiliary field, and derived field. Update observers.
void
pylith::feassemble::IntegratorInterface::initialize(const pylith::topology::Field& solution) {
    PYLITH_METHOD_BEGIN;
    PYLITH_JOURNAL_DEBUG("intialize(solution="<<solution.label()<<")");

    const bool isSubMesh = true;
    delete _interfaceMesh;_interfaceMesh = new pylith::topology::Mesh(isSubMesh);assert(_interfaceMesh);
    pylith::faults::TopologyOps::createFaultParallel(_interfaceMesh, solution.mesh(), _interfaceId,
                                                     _interfaceLabel.c_str());
    pylith::topology::MeshOps::checkTopology(*_interfaceMesh);

    // Optimize closure for coordinates.
    PetscDM dmFault = _interfaceMesh->dmMesh();assert(dmFault);
    pylith::topology::CoordsVisitor::optimizeClosure(dmFault);

    Integrator::initialize(solution);

    PYLITH_METHOD_END;
} // initialize


// ---------------------------------------------------------------------------------------------------------------------
// Compute RHS residual for G(t,s).
void
pylith::feassemble::IntegratorInterface::computeRHSResidual(pylith::topology::Field* residual,
                                                            const PylithReal t,
                                                            const PylithReal dt,
                                                            const pylith::topology::Field& solution) {
    PYLITH_METHOD_BEGIN;
    PYLITH_JOURNAL_DEBUG("computeRHSResidual(residual="<<residual<<", t="<<t<<", dt="<<dt<<", solution="<<solution.label()<<")");

    if ((0 == _kernelsRHSResidualPos.size()) && (0 == _kernelsRHSResidualNeg.size())) { PYLITH_METHOD_END;}

    _setKernelConstants(solution, dt);

    pylith::topology::Field solutionDot(solution.mesh()); // No dependence on time derivative of solution in RHS.
    solutionDot.label("solution_dot");

    const PylithInt faceDim = solution.mesh().dimension() - 1;

    const PylithInt labelValuePos = 100 + faceDim;
    _IntegratorInterface::computeResidual(residual, this, _kernelsRHSResidualPos, labelValuePos,
                                          t, dt, solution, solutionDot);

    const PylithInt labelValueNeg = -(100 + faceDim);
    _IntegratorInterface::computeResidual(residual, this, _kernelsRHSResidualPos, labelValueNeg,
                                          t, dt, solution, solutionDot);

    PYLITH_METHOD_END;
} // computeRHSResidual


// ---------------------------------------------------------------------------------------------------------------------
// Compute RHS Jacobian for G(t,s).
void
pylith::feassemble::IntegratorInterface::computeRHSJacobian(PetscMat jacobianMat,
                                                            PetscMat precondMat,
                                                            const PylithReal t,
                                                            const PylithReal dt,
                                                            const pylith::topology::Field& solution) {
    PYLITH_METHOD_BEGIN;
    PYLITH_JOURNAL_DEBUG("computeRHSJacobian(jacobianMat="<<jacobianMat<<", precondMat="<<precondMat<<", t="<<t<<", dt="<<dt<<", solution="<<solution.label()<<") empty method");

    if ((0 == _kernelsRHSJacobianPos.size()) && (0 == _kernelsRHSJacobianNeg.size())) { PYLITH_METHOD_END;}

    _setKernelConstants(solution, dt);

    pylith::topology::Field solutionDot(solution.mesh()); // No dependence on time derivative of solution in RHS.
    solutionDot.label("solution_dot");
    const PylithReal s_tshift = 0.0; // No dependence on time derivative of solution in RHS, so shift isn't applicable.

    const PylithInt faceDim = solution.mesh().dimension() - 1;

    const PylithInt labelValuePos = 100 + faceDim;
    _IntegratorInterface::computeJacobian(jacobianMat, precondMat, this, _kernelsRHSJacobianPos, labelValuePos,
                                          t, dt, s_tshift, solution, solutionDot);

    const PylithInt labelValueNeg = -(100 + faceDim);
    _IntegratorInterface::computeJacobian(jacobianMat, precondMat, this, _kernelsRHSJacobianPos, labelValueNeg,
                                          t, dt, s_tshift, solution, solutionDot);

    _needNewRHSJacobian = false;

    PYLITH_METHOD_END;
} // computeRHSJacobian


// ---------------------------------------------------------------------------------------------------------------------
// Compute LHS residual for F(t,s,\dot{s}).
void
pylith::feassemble::IntegratorInterface::computeLHSResidual(pylith::topology::Field* residual,
                                                            const PylithReal t,
                                                            const PylithReal dt,
                                                            const pylith::topology::Field& solution,
                                                            const pylith::topology::Field& solutionDot) {
    PYLITH_METHOD_BEGIN;
    PYLITH_JOURNAL_DEBUG("computeLHSResidual(residual="<<residual<<", t="<<t<<", dt="<<dt<<", solution="<<solution.label()<<")");

    if ((0 == _kernelsLHSResidualPos.size()) && (0 == _kernelsLHSResidualNeg.size())) { PYLITH_METHOD_END;}

    _setKernelConstants(solution, dt);

    const PylithInt faceDim = solution.mesh().dimension() - 1;

    const PylithInt labelValuePos = 100 + faceDim;
    _IntegratorInterface::computeResidual(residual, this, _kernelsLHSResidualPos, labelValuePos,
                                          t, dt, solution, solutionDot);

    const PylithInt labelValueNeg = -(100 + faceDim);
    _IntegratorInterface::computeResidual(residual, this, _kernelsLHSResidualPos, labelValueNeg,
                                          t, dt, solution, solutionDot);

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
    PYLITH_JOURNAL_DEBUG("computeLHSJacobian(jacobianMat="<<jacobianMat<<", precondMat="<<precondMat<<", t="<<t<<", dt="<<dt<<", solution="<<solution.label()<<", solutionDot="<<solutionDot.label()<<")");

    if ((0 == _kernelsLHSJacobianPos.size()) && (0 == _kernelsLHSJacobianNeg.size())) { PYLITH_METHOD_END;}

    _setKernelConstants(solution, dt);

    const PylithInt faceDim = solution.mesh().dimension() - 1;

    const PylithInt labelValuePos = 100 + faceDim;
    _IntegratorInterface::computeJacobian(jacobianMat, precondMat, this, _kernelsLHSJacobianPos, labelValuePos,
                                          t, dt, s_tshift, solution, solutionDot);

    const PylithInt labelValueNeg = -(100 + faceDim);
    _IntegratorInterface::computeJacobian(jacobianMat, precondMat, this, _kernelsLHSJacobianPos, labelValueNeg,
                                          t, dt, s_tshift, solution, solutionDot);

    _needNewLHSJacobian = false;

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
    PYLITH_JOURNAL_DEBUG("computeLHSJacobianLumpedInv(jacobianInv="<<jacobianInv<<", t="<<t<<", dt="<<dt<<", solution="<<solution.label()<<") empty method");

    // No implementation needed for interface.

    PYLITH_METHOD_END;
} // computeLHSJacobianLumpedInv


// ---------------------------------------------------------------------------------------------------------------------
// Compute residual.
void
pylith::feassemble::_IntegratorInterface::computeResidual(pylith::topology::Field* residual,
                                                          const pylith::feassemble::IntegratorInterface* integrator,
                                                          const std::vector<pylith::feassemble::IntegratorInterface::ResidualKernels>& kernels,
                                                          const PylithInt labelValue,
                                                          const PylithReal t,
                                                          const PylithReal dt,
                                                          const pylith::topology::Field& solution,
                                                          const pylith::topology::Field& solutionDot) {
    PYLITH_METHOD_BEGIN;
    // PYLITH_COMPONENT_DEBUG("computeRHSResidual(residual="<<residual<<", t="<<t<<", dt="<<dt<<",
    // solution="<<solution.label()<<")");

    assert(integrator);
    assert(residual);

    const pylith::topology::Field* auxiliaryField = integrator->getAuxiliaryField();assert(auxiliaryField);

    PetscErrorCode err;

    PetscDM dmSoln = solution.dmMesh();
    PetscDM dmAux = auxiliaryField->dmMesh();

    PetscDS prob = NULL;
    err = DMGetDS(dmSoln, &prob);PYLITH_CHECK_ERROR(err);

    // Get auxiliary data
    err = PetscObjectCompose((PetscObject) dmSoln, "dmAux", (PetscObject) dmAux);PYLITH_CHECK_ERROR(err);
    err = PetscObjectCompose((PetscObject) dmSoln, "A", (PetscObject) auxiliaryField->localVector());PYLITH_CHECK_ERROR(err);

    std::string interfaceLabelName = std::string(integrator->getSurfaceMarkerLabel()) + std::string("-interface");
    PetscDMLabel interfaceDMLabel = NULL;
    err = DMGetLabel(dmSoln, interfaceLabelName.c_str(), &interfaceDMLabel);PYLITH_CHECK_ERROR(err);assert(interfaceDMLabel);

#if 0 // DEBUGGING
    std::cout << "DOMAIN MESH" << std::endl;
    solution.mesh().view("::ascii_info_detail");

    std::cout << "FAULT MESH" << std::endl;
    auxiliaryField->mesh().view("::ascii_info_detail");

    err = DMLabelView(interfaceDMLabel, PETSC_VIEWER_STDOUT_SELF);
#endif // DEBUGGING

    // Compute the local residual
    assert(solution.localVector());
    assert(residual->localVector());
    for (size_t i = 0; i < kernels.size(); ++i) {
        const PetscInt i_field = solution.subfieldInfo(kernels[i].subfield.c_str()).index;
        err = PetscDSSetBdResidual(prob, i_field, kernels[i].r0, kernels[i].r1);PYLITH_CHECK_ERROR(err);
        err = DMPlexComputeBdResidualSingle(dmSoln, t, interfaceDMLabel, 1, &labelValue, i_field, solution.localVector(), solutionDot.localVector(), residual->localVector());PYLITH_CHECK_ERROR(err);
    } // for

    PYLITH_METHOD_END;
} // computeResidual


// ---------------------------------------------------------------------------------------------------------------------
// Compute Jacobian.
void
pylith::feassemble::_IntegratorInterface::computeJacobian(PetscMat jacobianMat,
                                                          PetscMat precondMat,
                                                          const pylith::feassemble::IntegratorInterface* integrator,
                                                          const std::vector<pylith::feassemble::IntegratorInterface::JacobianKernels>& kernels,
                                                          const PylithInt labelValue,
                                                          const PylithReal t,
                                                          const PylithReal dt,
                                                          const PylithReal s_tshift,
                                                          const pylith::topology::Field& solution,
                                                          const pylith::topology::Field& solutionDot) {
    PYLITH_METHOD_BEGIN;
    // PYLITH_JOURNAL_DEBUG("computeRHSJacobian(jacobianMat="<<jacobianMat<<", precondMat="<<precondMat<<", t="<<t<<",
    // dt="<<dt<<", solution="<<solution.label()<<")");

    assert(jacobianMat);
    assert(precondMat);

    const pylith::topology::Field* auxiliaryField = integrator->getAuxiliaryField();assert(auxiliaryField);

    PetscErrorCode err;

    PetscDM dmSoln = solution.dmMesh();
    PetscDM dmAux = auxiliaryField->dmMesh();

    PetscDS prob = NULL;
    err = DMGetDS(dmSoln, &prob);PYLITH_CHECK_ERROR(err);

    // Get auxiliary data
    err = PetscObjectCompose((PetscObject) dmSoln, "dmAux", (PetscObject) dmAux);PYLITH_CHECK_ERROR(err);
    err = PetscObjectCompose((PetscObject) dmSoln, "A", (PetscObject) auxiliaryField->localVector());PYLITH_CHECK_ERROR(err);

    std::string interfaceLabelName = std::string(integrator->getSurfaceMarkerLabel()) + std::string("-interface");
    PetscDMLabel interfaceDMLabel = NULL;
    err = DMGetLabel(dmSoln, interfaceLabelName.c_str(), &interfaceDMLabel);PYLITH_CHECK_ERROR(err);assert(interfaceDMLabel);

    // Compute the local Jacobian
    assert(solution.localVector());

    for (size_t i = 0; i < kernels.size(); ++i) {
        const PetscInt i_fieldTrial = solution.subfieldInfo(kernels[i].subfieldTrial.c_str()).index;
        const PetscInt i_fieldBasis = solution.subfieldInfo(kernels[i].subfieldBasis.c_str()).index;
        err = PetscDSSetBdJacobian(prob, i_fieldTrial, i_fieldBasis, kernels[i].j0, kernels[i].j1, kernels[i].j2, kernels[i].j3);PYLITH_CHECK_ERROR(err);

        err = DMPlexComputeBdJacobianSingle(dmSoln, t, interfaceDMLabel, 1, &labelValue, i_fieldTrial, solution.localVector(),
                                            solutionDot.localVector(), s_tshift, jacobianMat, precondMat);PYLITH_CHECK_ERROR(err);
    } // for

    PYLITH_METHOD_END;
} // computeJacobian


// End of file
