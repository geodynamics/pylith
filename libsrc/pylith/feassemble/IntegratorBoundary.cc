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

#include "IntegratorBoundary.hh" // implementation of object methods

#include "pylith/topology/Mesh.hh" // USES Mesh
#include "pylith/topology/MeshOps.hh" // USES createLowerDimMesh()
#include "pylith/topology/Field.hh" // USES Field
#include "pylith/problems/Physics.hh" // USES Physics
#include "pylith/topology/CoordsVisitor.hh" // USES CoordsVisitor::optimizeClosure()

#include "spatialdata/spatialdb/GravityField.hh" // HASA GravityField

#include "pylith/utils/journals.hh" // USES PYLITH_JOURNAL_*

#include "petscds.h" // USES PetscDS

#include <cassert> // USES assert()
#include <stdexcept> // USES std::runtime_error

// ---------------------------------------------------------------------------------------------------------------------
// Local "private" functions.
namespace pylith {
    namespace feassemble {
        class _IntegratorBoundary {
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
                                 const pylith::feassemble::IntegratorBoundary* integrator,
                                 const std::vector<pylith::feassemble::IntegratorBoundary::ResidualKernels>& kernels,
                                 const PylithReal t,
                                 const PylithReal dt,
                                 const pylith::topology::Field& solution,
                                 const pylith::topology::Field& solutionDot);

            static const char* genericComponent;
        }; // _IntegratorBoundary
        const char* _IntegratorBoundary::genericComponent = "integratorboundary";

    } // feassemble
} // pylith

// ---------------------------------------------------------------------------------------------------------------------
// Default constructor.
pylith::feassemble::IntegratorBoundary::IntegratorBoundary(pylith::problems::Physics* const physics) :
    Integrator(physics),
    _boundaryMesh(NULL),
    _boundaryLabel("") {
    GenericComponent::setName(_IntegratorBoundary::genericComponent);
} // constructor


// ---------------------------------------------------------------------------------------------------------------------
// Destructor.
pylith::feassemble::IntegratorBoundary::~IntegratorBoundary(void) {
    deallocate();
} // destructor


// ---------------------------------------------------------------------------------------------------------------------
// Deallocate PETSc and local data structures.
void
pylith::feassemble::IntegratorBoundary::deallocate(void) {
    PYLITH_METHOD_BEGIN;

    Integrator::deallocate();

    delete _boundaryMesh;_boundaryMesh = NULL;

    PYLITH_METHOD_END;
} // deallocate


// ---------------------------------------------------------------------------------------------------------------------
// Set label marking boundary associated with boundary condition surface.
void
pylith::feassemble::IntegratorBoundary::setMarkerLabel(const char* value) {
    PYLITH_JOURNAL_DEBUG("setMarkerLabel(value="<<value<<")");

    if (strlen(value) == 0) {
        throw std::runtime_error("Empty string given for boundary condition integrator label.");
    } // if

    _boundaryLabel = value;
} // setMarkerLabel


// ---------------------------------------------------------------------------------------------------------------------
// Get label marking boundary associated with boundary condition surface.
const char*
pylith::feassemble::IntegratorBoundary::getMarkerLabel(void) const {
    return _boundaryLabel.c_str();
} // getMarkerLabel


// ---------------------------------------------------------------------------------------------------------------------
// Get mesh associated with integration domain.
const pylith::topology::Mesh&
pylith::feassemble::IntegratorBoundary::getPhysicsDomainMesh(void) const {
    assert(_boundaryMesh);
    return *_boundaryMesh;
} // domainMesh


// ---------------------------------------------------------------------------------------------------------------------
void
pylith::feassemble::IntegratorBoundary::setKernelsRHSResidual(const std::vector<ResidualKernels>& kernels) {
    PYLITH_METHOD_BEGIN;
    PYLITH_JOURNAL_DEBUG("setKernelsRHSResidual(# kernels="<<kernels.size()<<")");

    _kernelsRHSResidual = kernels;

    PYLITH_METHOD_END;
} // setKernelsRHSResidual


// ---------------------------------------------------------------------------------------------------------------------
void
pylith::feassemble::IntegratorBoundary::setKernelsLHSResidual(const std::vector<ResidualKernels>& kernels) {
    PYLITH_METHOD_BEGIN;
    PYLITH_JOURNAL_DEBUG("setKernelsLHSResidual(# kernels="<<kernels.size()<<")");

    _kernelsLHSResidual = kernels;

    PYLITH_METHOD_END;
} // setKernelsLHSResidual


// ---------------------------------------------------------------------------------------------------------------------
// Initialize integration domain, auxiliary field, and derived field. Update observers.
void
pylith::feassemble::IntegratorBoundary::initialize(const pylith::topology::Field& solution) {
    PYLITH_METHOD_BEGIN;
    PYLITH_JOURNAL_DEBUG("intialize(solution="<<solution.getLabel()<<")");

    delete _boundaryMesh;
    _boundaryMesh = pylith::topology::MeshOps::createLowerDimMesh(solution.mesh(), _boundaryLabel.c_str());
    assert(_boundaryMesh);
    pylith::topology::CoordsVisitor::optimizeClosure(_boundaryMesh->dmMesh());

    Integrator::initialize(solution);

    PYLITH_METHOD_END;
} // initialize


// ---------------------------------------------------------------------------------------------------------------------
// Update auxiliary field values to current time.
void
pylith::feassemble::IntegratorBoundary::updateState(const double t) {
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
pylith::feassemble::IntegratorBoundary::computeRHSResidual(pylith::topology::Field* residual,
                                                           const PylithReal t,
                                                           const PylithReal dt,
                                                           const pylith::topology::Field& solution) {
    PYLITH_METHOD_BEGIN;
    PYLITH_JOURNAL_DEBUG("computeRHSResidual(residual="<<residual<<", t="<<t<<", dt="<<dt<<", solution="<<solution.getLabel()<<")");

    if (0 == _kernelsRHSResidual.size()) { PYLITH_METHOD_END;}

    _setKernelConstants(solution, dt);

    pylith::topology::Field solutionDot(solution.mesh()); // No dependence on time derivative of solution in RHS.
    solutionDot.setLabel("solution_dot");
    _IntegratorBoundary::computeResidual(residual, this, _kernelsRHSResidual, t, dt, solution, solutionDot);

    PYLITH_METHOD_END;
} // computeRHSResidual


// ---------------------------------------------------------------------------------------------------------------------
// Compute RHS Jacobian for G(t,s).
void
pylith::feassemble::IntegratorBoundary::computeRHSJacobian(PetscMat jacobianMat,
                                                           PetscMat precondMat,
                                                           const PylithReal t,
                                                           const PylithReal dt,
                                                           const pylith::topology::Field& solution) {
    PYLITH_METHOD_BEGIN;
    PYLITH_JOURNAL_DEBUG("computeRHSJacobian(jacobianMat="<<jacobianMat<<", precondMat="<<precondMat<<", t="<<t<<", dt="<<dt<<", solution="<<solution.getLabel()<<") empty method");

    _needNewRHSJacobian = false;
    // No implementation needed for boundary.

    PYLITH_METHOD_END;
} // computeRHSJacobian


// ---------------------------------------------------------------------------------------------------------------------
// Compute LHS residual for F(t,s,\dot{s}).
void
pylith::feassemble::IntegratorBoundary::computeLHSResidual(pylith::topology::Field* residual,
                                                           const PylithReal t,
                                                           const PylithReal dt,
                                                           const pylith::topology::Field& solution,
                                                           const pylith::topology::Field& solutionDot) {
    PYLITH_METHOD_BEGIN;
    PYLITH_JOURNAL_DEBUG("computeLHSResidual(residual="<<residual<<", t="<<t<<", dt="<<dt<<", solution="<<solution.getLabel()<<")");

    if (0 == _kernelsLHSResidual.size()) { PYLITH_METHOD_END;}

    _setKernelConstants(solution, dt);

    _IntegratorBoundary::computeResidual(residual, this, _kernelsLHSResidual, t, dt, solution, solutionDot);
} // computeLHSResidual


// ---------------------------------------------------------------------------------------------------------------------
// Compute LHS Jacobian for F(t,s,\dot{s}).
void
pylith::feassemble::IntegratorBoundary::computeLHSJacobian(PetscMat jacobianMat,
                                                           PetscMat precondMat,
                                                           const PylithReal t,
                                                           const PylithReal dt,
                                                           const PylithReal s_tshift,
                                                           const pylith::topology::Field& solution,
                                                           const pylith::topology::Field& solutionDot) {
    PYLITH_METHOD_BEGIN;
    PYLITH_JOURNAL_DEBUG("computeLHSJacobian(jacobianMat="<<jacobianMat<<", precondMat="<<precondMat<<", t="<<t<<", dt="<<dt<<", solution="<<solution.getLabel()<<", solutionDot="<<solutionDot.getLabel()<<") empty method");

    _needNewLHSJacobian = false;
    // No implementation needed for boundary.

    PYLITH_METHOD_END;
} // computeLHSJacobian


// ---------------------------------------------------------------------------------------------------------------------
// Compute LHS Jacobian for F(t,s,\dot{s}).
void
pylith::feassemble::IntegratorBoundary::computeLHSJacobianLumpedInv(pylith::topology::Field* jacobianInv,
                                                                    const PylithReal t,
                                                                    const PylithReal dt,
                                                                    const PylithReal s_tshift,
                                                                    const pylith::topology::Field& solution) {
    PYLITH_METHOD_BEGIN;
    PYLITH_JOURNAL_DEBUG("computeLHSJacobianLumpedInv(jacobianInv="<<jacobianInv<<", t="<<t<<", dt="<<dt<<", solution="<<solution.getLabel()<<") empty method");

    _needNewLHSJacobianLumped = false;
    // No implementation needed for boundary.

    PYLITH_METHOD_END;
} // computeLHSJacobianLumpedInv


// ---------------------------------------------------------------------------------------------------------------------
// Compute residual.
void
pylith::feassemble::_IntegratorBoundary::computeResidual(pylith::topology::Field* residual,
                                                         const pylith::feassemble::IntegratorBoundary* integrator,
                                                         const std::vector<pylith::feassemble::IntegratorBoundary::ResidualKernels>& kernels,
                                                         const PylithReal t,
                                                         const PylithReal dt,
                                                         const pylith::topology::Field& solution,
                                                         const pylith::topology::Field& solutionDot) {
    PYLITH_METHOD_BEGIN;
    pythia::journal::debug_t debug(_IntegratorBoundary::genericComponent);
    debug << pythia::journal::at(__HERE__)
          << "_IntegratorBoundary::computeRHSResidual(residual="<<residual<<", integrator="<<integrator
          <<", # kernels="<<kernels.size()<<", t="<<t<<", dt="<<dt<<", solution="<<solution.getLabel()
          <<", solutionDot="<<solutionDot.getLabel()<<")"
          << pythia::journal::endl;

    assert(integrator);
    assert(residual);

    const pylith::topology::Field* auxiliaryField = integrator->getAuxiliaryField();assert(auxiliaryField);

    PetscErrorCode err;

    // :KLUDGE: Potentially we may have multiple PetscDS objects. This assumes that the first one (with a NULL label) is
    // the correct one.
    PetscDS prob = NULL;
    PetscDM dmSoln = solution.dmMesh();assert(dmSoln);
    err = DMGetDS(dmSoln, &prob);PYLITH_CHECK_ERROR(err);assert(prob);

    // Get auxiliary data
    PetscDM dmAux = auxiliaryField->dmMesh();assert(dmAux);
    err = PetscObjectCompose((PetscObject) dmSoln, "dmAux", (PetscObject) dmAux);PYLITH_CHECK_ERROR(err);
    err = PetscObjectCompose((PetscObject) dmSoln, "A", (PetscObject) auxiliaryField->localVector());PYLITH_CHECK_ERROR(err);

    // Compute the local residual
    assert(solution.localVector());
    assert(residual->localVector());
    PetscDMLabel dmLabel = NULL;
    err = DMGetLabel(dmSoln, integrator->getMarkerLabel(), &dmLabel);PYLITH_CHECK_ERROR(err);
    const int labelId = 1;

    // solution.mesh().view(":mesh.txt:ascii_info_detail"); // :DEBUG:

    for (size_t i = 0; i < kernels.size(); ++i) {
        const PetscInt i_field = solution.subfieldInfo(kernels[i].subfield.c_str()).index;
        err = PetscDSSetBdResidual(prob, i_field, kernels[i].r0, kernels[i].r1);PYLITH_CHECK_ERROR(err);
        err = DMPlexComputeBdResidualSingle(dmSoln, t, dmLabel, 1, &labelId, i_field, solution.localVector(), solutionDot.localVector(), residual->localVector());PYLITH_CHECK_ERROR(err);
    } // for

    PYLITH_METHOD_END;
} // _computeResidual


// End of file
