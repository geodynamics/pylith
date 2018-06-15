// -*- C++ -*-
//
// ----------------------------------------------------------------------
//
// Brad T. Aagaard, U.S. Geological Survey
// Charles A. Williams, GNS Science
// Matthew G. Knepley, University of Chicago
//
// This code was developed as part of the Computational Infrastructure
// for Geodynamics (http://geodynamics.org).
//
// Copyright (c) 2010-2016 University of California, Davis
//
// See COPYING for license information.
//
// ----------------------------------------------------------------------
//

#include <portinfo>

#include "IntegratorBoundary.hh" // implementation of object methods

#include "pylith/feassemble/AuxiliaryFactory.hh" // USES AuxiliaryFactory
#include "pylith/topology/Mesh.hh" // USES Mesh
#include "pylith/topology/Field.hh" // USES Field
#include "pylith/topology/CoordsVisitor.hh" // USES CoordsVisitor

#include "pylith/utils/journals.hh" // USES PYLITH_COMPONENT_*
#include "pylith/utils/error.hh" // USES PYLITH_METHOD_BEGIN/END

#include <strings.h> // USES strcasecmp()
#include <cassert> // USES assert()
#include <stdexcept> // USES std::runtime_error
#include <sstream> \
    // USES std::ostringstream

// ----------------------------------------------------------------------
// Local "private" functions.
namespace pylith {
    namespace feassemble {
        static
        void _computeResidual(pylith::topology::Field* residual,
                              const pylith::feassemble::IntegratorBoundary* integrator,
                              const PylithReal t,
                              const PylithReal dt,
                              const pylith::topology::Field& solution,
                              const pylith::topology::Field& solutionDot);
    } // feassemble
} // pylith

// ----------------------------------------------------------------------
// Default constructor.
pylith::feassemble::IntegratorBoundary::IntegratorBoundary(void) :
    _boundaryMesh(NULL)
{ // constructor
    _refDir1[0] = 0.0;
    _refDir1[1] = 0.0;
    _refDir1[2] = 1.0;

    _refDir2[0] = 0.0;
    _refDir2[1] = 1.0;
    _refDir2[2] = 0.0;
} // constructor


// ----------------------------------------------------------------------
// Destructor.
pylith::feassemble::IntegratorBoundary::~IntegratorBoundary(void) {
    deallocate();
} // destructor


// ----------------------------------------------------------------------
// Deallocate PETSc and local data structures.
void
pylith::feassemble::IntegratorBoundary::deallocate(void) {
    PYLITH_METHOD_BEGIN;

    IntegratorPointwise::deallocate();

    delete _boundaryMesh; _boundaryMesh = NULL;

    PYLITH_METHOD_END;
} // deallocate


// ----------------------------------------------------------------------
// Set mesh label associated with boundary condition surface.
void
pylith::feassemble::IntegratorBoundary::label(const char* value) {
    if (strlen(value) == 0) {
        throw std::runtime_error("Empty string given for boundary condition label.");
    } // if

    _label = value;
} // label


// ----------------------------------------------------------------------
// Get mesh label associated with boundary condition surface.
const char*
pylith::feassemble::IntegratorBoundary::label(void) const {
    return _label.c_str();
} // label


// ----------------------------------------------------------------------
// Set name of field in solution to constrain.
void
pylith::feassemble::IntegratorBoundary::field(const char* value) {
    PYLITH_METHOD_BEGIN;

    if (strlen(value) == 0) {
        throw std::runtime_error("Empty string given for name of solution field for boundary condition.");
    } // if
    _field = value;

    PYLITH_METHOD_END;
}  // field


// ----------------------------------------------------------------------
// Get name of field in solution to constrain.
const char*
pylith::feassemble::IntegratorBoundary::field(void) const {
    journal::debug_t debug("boundarycondition");
    debug << journal::at(__HERE__)
          << "BoundaryCondition::field()" << journal::endl;

    return _field.c_str();
} // field


// ----------------------------------------------------------------------
// Set first choice for reference direction to discriminate among tangential directions in 3-D.
void
pylith::feassemble::IntegratorBoundary::refDir1(const PylithReal vec[3]) {
    // Set reference direction, insuring it is a unit vector.
    const PylithReal mag = sqrt(vec[0]*vec[0] + vec[1]*vec[1] + vec[2]*vec[2]);
    if (mag < 1.0e-6) {
        std::ostringstream msg;
        msg << "Magnitude of reference direction 1 ("<<vec[0]<<", "<<vec[1]<<", "<<vec[2]
            <<") for boundary condition '" << label() << "' is negligible. Use a unit vector.";
        throw std::runtime_error(msg.str());
    } // if
    for (int i = 0; i < 3; ++i) {
        _refDir1[i] = vec[i] / mag;
    } // for
} // refDir1


// ----------------------------------------------------------------------
// Set second choice for reference direction to discriminate among tangential directions in 3-D.
void
pylith::feassemble::IntegratorBoundary::refDir2(const PylithReal vec[3]) {
    // Set reference direction, insuring it is a unit vector.
    const PylithReal mag = sqrt(vec[0]*vec[0] + vec[1]*vec[1] + vec[2]*vec[2]);
    if (mag < 1.0e-6) {
        std::ostringstream msg;
        msg << "Magnitude of reference direction 2 ("<<vec[0]<<", "<<vec[1]<<", "<<vec[2]
            <<") for boundary condition '" << label() << "' is negligible. Use a unit vector.";
        throw std::runtime_error(msg.str());
    } // if
    for (int i = 0; i < 3; ++i) {
        _refDir1[i] = vec[i] / mag;
    } // for
} // refDir2


// ----------------------------------------------------------------------
// Get mesh associated with integrator domain.
const pylith::topology::Mesh&
pylith::feassemble::IntegratorBoundary::domainMesh(void) const {
    assert(_boundaryMesh);
    return *_boundaryMesh;
} // domainMesh


// ----------------------------------------------------------------------
// Verify configuration is acceptable.
void
pylith::feassemble::IntegratorBoundary::verifyConfiguration(const pylith::topology::Field& solution) const {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("verifyConfiguration(solution="<<solution.label()<<")");

    if (!solution.hasSubfield(_field.c_str())) {
        std::ostringstream msg;
        msg << "Cannot apply Neumann boundary condition to field '"<< _field
            << "'; field is not in solution.";
        throw std::runtime_error(msg.str());
    } // if

    PYLITH_METHOD_END;
} // verifyConfiguration


// ----------------------------------------------------------------------
// Initialize boundary condition.
void
pylith::feassemble::IntegratorBoundary::initialize(const pylith::topology::Field& solution) {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("initialize(solution="<<solution.label()<<")");

    _boundaryMesh = new pylith::topology::Mesh(solution.mesh(), _label.c_str()); assert(_boundaryMesh);
    PetscDM dmBoundary = _boundaryMesh->dmMesh(); assert(dmBoundary);
    pylith::topology::CoordsVisitor::optimizeClosure(dmBoundary);

    delete _auxField; _auxField = new pylith::topology::Field(*_boundaryMesh); assert(_auxField);
    _auxField->label("auxiliary subfields");
    _auxFieldSetup(solution);
    _auxField->subfieldsSetup();
    _auxField->allocate();
    _auxField->zeroLocal();

    assert(_normalizer);
    pylith::feassemble::AuxiliaryFactory* factory = _auxFactory(); assert(factory);
    factory->initializeSubfields();

    //_auxField->view("AUXILIARY FIELD"); // :DEBUG:

    delete _derivedField; _derivedField = NULL;

    const bool infoOnly = true;
    notifyObservers(0.0, 0, solution, infoOnly);

    PYLITH_METHOD_END;
} // initialize

// ----------------------------------------------------------------------
// Compute RHS residual for G(t,s).
void
pylith::feassemble::IntegratorBoundary::computeRHSResidual(pylith::topology::Field* residual,
                                                           const PylithReal t,
                                                           const PylithReal dt,
                                                           const pylith::topology::Field& solution) {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("computeRHSResidual(residual="<<residual<<", t="<<t<<", dt="<<dt<<", solution="<<solution.label()<<")");

    if (!_hasFEKernels(KERNELS_RHS_RESIDUAL)) { PYLITH_METHOD_END; }

    pylith::topology::Field solutionDot(solution.mesh()); // No dependence on time derivative of solution in RHS.
    solutionDot.label("solution_dot");

    _setFEKernels(solution, KERNELS_RHS_RESIDUAL);
    _setFEConstants(solution, dt);
    _computeResidual(residual, this, t, dt, solution, solutionDot);

    PYLITH_METHOD_END;
} // computeRHSResidual


// ----------------------------------------------------------------------
// Compute RHS residual for G(t,s).
void
pylith::feassemble::IntegratorBoundary::computeLHSResidual(pylith::topology::Field* residual,
                                                           const PylithReal t,
                                                           const PylithReal dt,
                                                           const pylith::topology::Field& solution,
                                                           const pylith::topology::Field& solutionDot) {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("computeLHSResidual(residual="<<residual<<", t="<<t<<", dt="<<dt<<", solution="<<solution.label()<<", solutionDot="<<solutionDot.label()<<")");

    if (!_hasFEKernels(KERNELS_LHS_RESIDUAL)) { PYLITH_METHOD_END; }

    _setFEKernels(solution, KERNELS_LHS_RESIDUAL);
    _setFEConstants(solution, dt);
    _computeResidual(residual, this, t, dt, solution, solutionDot);

    PYLITH_METHOD_END;
} // computeRHSResidual


// ----------------------------------------------------------------------
// Compute RHS Jacobian and preconditioner for G(t,s).
void
pylith::feassemble::IntegratorBoundary::computeRHSJacobian(PetscMat jacobianMat,
                                                           PetscMat preconMat,
                                                           const PylithReal t,
                                                           const PylithReal dt,
                                                           const pylith::topology::Field& solution) {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("computeRHSJacobian(jacobianMat="<<jacobianMat<<", preconMat="<<preconMat<<", t="<<t<<", dt="<<dt<<", solution="<<solution.label()<<")");

    if (!_hasFEKernels(KERNELS_RHS_JACOBIAN)) { PYLITH_METHOD_END; }

    _setFEKernels(solution, KERNELS_RHS_JACOBIAN);
    _setFEConstants(solution, dt);
    PYLITH_COMPONENT_ERROR(":TODO: Implement IntegratorBoundary::_computeJacobian().");
    throw std::logic_error(":TODO: Implement IntegratorBoundary::_computeJacobian().");

    PYLITH_METHOD_END;
} // computeRHSJacobian


// ----------------------------------------------------------------------
// Compute LHS Jacobian and preconditioner for F(t,s,\dot{s}) with implicit time-stepping.
void
pylith::feassemble::IntegratorBoundary::computeLHSJacobianImplicit(PetscMat jacobianMat,
                                                                   PetscMat precondMat,
                                                                   const PylithReal t,
                                                                   const PylithReal dt,
                                                                   const PylithReal s_tshift,
                                                                   const pylith::topology::Field& solution,
                                                                   const pylith::topology::Field& solutionDot) {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("computeLHSJacobianImplicit(jacobianMat="<<jacobianMat<<", precondMat="<<precondMat<<", t="<<t<<", dt="<<dt<<", solution="<<solution.label()<<", solutionDot="<<solutionDot.label()<<")");

    assert(jacobianMat);
    assert(precondMat);
    assert(s_tshift > 0);

    if (!_hasFEKernels(KERNELS_LHS_JACOBIAN)) { PYLITH_METHOD_END; }

    _setFEKernels(solution, KERNELS_LHS_JACOBIAN);
    _setFEConstants(solution, dt);
    PYLITH_COMPONENT_ERROR(":TODO: Implement IntegratorBoundary::_computeJacobian().");
    throw std::logic_error(":TODO: Implement IntegratorBoundary::_computeJacobian().");

    PYLITH_METHOD_END;
} // computeLHSJacobianImplicit


// ----------------------------------------------------------------------
// Compute inverse of lumped LHS Jacobian for F(t,s,\dot{s}) with explicit time-stepping.
void
pylith::feassemble::IntegratorBoundary::computeLHSJacobianLumpedInv(pylith::topology::Field* jacobianInv,
                                                                    const PylithReal t,
                                                                    const PylithReal dt,
                                                                    const PylithReal s_tshift,
                                                                    const pylith::topology::Field& solution) {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("computeLHSJacobianLumpedInv(jacobianMat="<<jacobianInv<<", t="<<t<<", dt="<<dt<<", s_tshift="<<s_tshift<<", solution="<<solution.label()<<")");

    assert(jacobianInv);
    assert(s_tshift > 0);

    if (!_hasFEKernels(KERNELS_LHS_JACOBIAN_LUMPEDINV)) { PYLITH_METHOD_END; }

    _setFEKernels(solution, KERNELS_LHS_JACOBIAN_LUMPEDINV);
    _setFEConstants(solution, dt);
    PYLITH_COMPONENT_ERROR(":TODO: Implement IntegratorBoundary::_computeJacobian().");
    throw std::logic_error(":TODO: Implement IntegratorBoundary::_computeJacobian().");

    PYLITH_METHOD_END;
} // computeLHSJacobianLumpedInv


// ----------------------------------------------------------------------
// Set constants used in finite-element integrations.
void
pylith::feassemble::IntegratorBoundary::_setFEConstants(const pylith::topology::Field& solution,
                                                        const PylithReal dt) const {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("_setFEConstants(solution="<<solution.label()<<", dt="<<dt<<")");

    const PetscInt numConstants = 6;
    PetscScalar constants[6];
    int index = 0;
    for (int i = 0; i < 3; ++i, index++) {
        constants[index] = _refDir1[i];
    } // for
    for (int i = 0; i < 3; ++i, index++) {
        constants[index] = _refDir2[i];
    } // for

    // Put up direction into constants.
    const PetscDM dmSoln = solution.dmMesh(); assert(dmSoln);
    PetscDS prob = NULL;
    PetscErrorCode err = DMGetDS(dmSoln, &prob); PYLITH_CHECK_ERROR(err); assert(prob);
    err = PetscDSSetConstants(prob, numConstants, constants); PYLITH_CHECK_ERROR(err);

    PYLITH_METHOD_END;
} // _setFEConstants


// ----------------------------------------------------------------------
// Compute residual.
void
pylith::feassemble::_computeResidual(pylith::topology::Field* residual,
                                     const pylith::feassemble::IntegratorBoundary* integrator,
                                     const PylithReal t,
                                     const PylithReal dt,
                                     const pylith::topology::Field& solution,
                                     const pylith::topology::Field& solutionDot) {
    PYLITH_METHOD_BEGIN;
    //PYLITH_COMPONENT_DEBUG("computeRHSResidual(residual="<<residual<<", t="<<t<<", dt="<<dt<<", solution="<<solution.label()<<")");

    assert(integrator);
    assert(residual);

    const pylith::topology::Field* auxField = integrator->auxField();assert(auxField);

    PetscDS prob = NULL;
    PetscErrorCode err;

    PetscDM dmSoln = solution.dmMesh();
    PetscDM dmAux = auxField->dmMesh();
    PetscDMLabel dmLabel = NULL;

    // Pointwise function have been set in DS
    err = DMGetDS(dmSoln, &prob); PYLITH_CHECK_ERROR(err);

    // Get auxiliary data
    err = PetscObjectCompose((PetscObject) dmSoln, "dmAux", (PetscObject) dmAux); PYLITH_CHECK_ERROR(err);
    err = PetscObjectCompose((PetscObject) dmSoln, "A", (PetscObject) auxField->localVector()); PYLITH_CHECK_ERROR(err);

    // Compute the local residual
    assert(solution.localVector());
    assert(residual->localVector());
    err = DMGetLabel(dmSoln, integrator->label(), &dmLabel); PYLITH_CHECK_ERROR(err);
    const int labelId = 1;
    const topology::Field::SubfieldInfo& info = solution.subfieldInfo(integrator->field());

    //solution.mesh().view(":mesh.txt:ascii_info_detail"); // :DEBUG:

    //PYLITH_COMPONENT_DEBUG("DMPlexComputeBdResidualSingle() for boundary '"<<label()<<"')");
    err = DMPlexComputeBdResidualSingle(dmSoln, t, dmLabel, 1, &labelId, info.index, solution.localVector(), solutionDot.localVector(), residual->localVector()); PYLITH_CHECK_ERROR(err);

    PYLITH_METHOD_END;
} // _computeResidual


// End of file
