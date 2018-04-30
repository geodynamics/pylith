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

#include "Neumann.hh" // implementation of object methods

#include "pylith/topology/Mesh.hh" // USES Mesh
#include "pylith/topology/Field.hh" // USES Field
#include "pylith/feassemble/AuxiliaryFactory.hh" // USES AuxiliaryFactory
#include "pylith/topology/VisitorMesh.hh" // USES VecVisitorMesh
#include "pylith/topology/FieldQuery.hh" // HOLDSA FieldQuery
#include "pylith/topology/CoordsVisitor.hh" // USES CoordsVisitor
#include "pylith/meshio/OutputManager.hh" // USES OutputManager
#include "spatialdata/geocoords/CoordSys.hh" // USES CoordSys
#include "spatialdata/units/Nondimensional.hh" // USES Nondimensional

#include "pylith/utils/journals.hh" // USES PYLITH_COMPONENT_*

#include <strings.h> // USES strcasecmp()
#include <cassert> // USES assert()
#include <stdexcept> // USES std::runtime_error
#include <sstream> \
    // USES std::ostringstream

// ----------------------------------------------------------------------
// Default constructor.
pylith::bc::Neumann::Neumann(void) :
    _boundaryMesh(NULL),
    _scaleName("pressure")
{ // constructor
    _refDir1[0] = 0.0;
    _refDir1[1] = 0.0;
    _refDir1[2] = 1.0;

    _refDir2[0] = 0.0;
    _refDir2[1] = 1.0;
    _refDir2[2] = 0.0;

    _description.label = "unknown";
    _description.vectorFieldType = pylith::topology::FieldBase::OTHER;
    _description.numComponents = 0;
    _description.scale = 1.0;
    _description.validator = NULL;

} // constructor

// ----------------------------------------------------------------------
// Destructor.
pylith::bc::Neumann::~Neumann(void) {
    deallocate();
} // destructor

// ----------------------------------------------------------------------
// Deallocate PETSc and local data structures.
void
pylith::bc::Neumann::deallocate(void) {
    PYLITH_METHOD_BEGIN;

    IntegratorPointwise::deallocate();

    delete _boundaryMesh; _boundaryMesh = NULL;

    PYLITH_METHOD_END;
} // deallocate

// ----------------------------------------------------------------------
// Name of scale associated with Neumann boundary condition (e.g., pressure for elasticity).
void
pylith::bc::Neumann::scaleName(const char* value) {
    if (( value == std::string("length")) ||
        ( value == std::string("time")) ||
        ( value == std::string("pressure")) ||
        ( value == std::string("density")) ||
        ( value == std::string("pressure")) ) {
        _scaleName = value;
    } else {
        std::ostringstream msg;
        msg << "Unknown name of scale ("<<value<<") for Neumann boundary condition '" << label() << "'.";
        throw std::runtime_error(msg.str());
    } // if
} // scaleName

// ----------------------------------------------------------------------
// Set first choice for reference direction to discriminate among tangential directions in 3-D.
void
pylith::bc::Neumann::refDir1(const double vec[3]) {
    // Set reference direction, insuring it is a unit vector.
    const PylithReal mag = sqrt(vec[0]*vec[0] + vec[1]*vec[1] + vec[2]*vec[2]);
    if (mag < 1.0e-6) {
        std::ostringstream msg;
        msg << "Magnitude of reference direction 1 ("<<vec[0]<<", "<<vec[1]<<", "<<vec[2]
            <<") for Neumann boundary condition '" << label() << "' is negligible. Use a unit vector.";
        throw std::runtime_error(msg.str());
    } // if
    for (int i = 0; i < 3; ++i) {
        _refDir1[i] = vec[i] / mag;
    } // for
} // refDir1

// ----------------------------------------------------------------------
// Set second choice for reference direction to discriminate among tangential directions in 3-D.
void
pylith::bc::Neumann::refDir2(const double vec[3]) {
    // Set reference direction, insuring it is a unit vector.
    const PylithReal mag = sqrt(vec[0]*vec[0] + vec[1]*vec[1] + vec[2]*vec[2]);
    if (mag < 1.0e-6) {
        std::ostringstream msg;
        msg << "Magnitude of reference direction 2 ("<<vec[0]<<", "<<vec[1]<<", "<<vec[2]
            <<") for Neumann boundary condition '" << label() << "' is negligible. Use a unit vector.";
        throw std::runtime_error(msg.str());
    } // if
    for (int i = 0; i < 3; ++i) {
        _refDir1[i] = vec[i] / mag;
    } // for
} // refDir2

// ----------------------------------------------------------------------
// Verify configuration is acceptable.
void
pylith::bc::Neumann::verifyConfiguration(const pylith::topology::Field& solution) const {
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
pylith::bc::Neumann::initialize(const pylith::topology::Field& solution) {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("initialize(solution="<<solution.label()<<")");

    const topology::Field::SubfieldInfo& info = solution.subfieldInfo(_field.c_str());
    _description = info.description;

    _setFEKernelsRHSResidual(solution);

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

    _auxField->view("AUXILIARY FIELD"); // :DEBUG:
    if (_output) {
        _output->writeInfo(*_auxField);
    } // if

#if 0 // Use low-level function to set kernels
    const PetscDM dmSoln = solution.dmMesh(); assert(dmSoln);
    PetscDS prob = NULL;
    PetscErrorCode err = DMGetDS(dmSoln, &prob); PYLITH_CHECK_ERROR(err);

    // :TODO: @brad @matt Expect this to need updating once we associate point functions with a domain.
    void* context = NULL;
    const int labelId = 1;
    const topology::Field::SubfieldInfo& info = solution.subfieldInfo(_field.c_str());
    err = PetscDSAddBoundary(prob, DM_BC_NATURAL, label(), label(), info.index, 0, NULL,
                             NULL, 1, &labelId, context); PYLITH_CHECK_ERROR(err);
#endif

    PYLITH_METHOD_END;
} // initialize

// ----------------------------------------------------------------------
// Compute RHS residual for G(t,s).
void
pylith::bc::Neumann::computeRHSResidual(pylith::topology::Field* residual,
                                        const PylithReal t,
                                        const PylithReal dt,
                                        const pylith::topology::Field& solution) {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("computeRHSResidual(residual="<<residual<<", t="<<t<<", dt="<<dt<<", solution="<<solution.label()<<")");

    _setFEKernelsRHSResidual(solution);
    _setFEConstants(solution, dt);

    pylith::topology::Field solutionDot(solution.mesh()); // No dependence on time derivative of solution in RHS.
    solutionDot.label("solution_dot");

    assert(residual);
    assert(_auxField);

    PetscDS prob = NULL;
    PetscErrorCode err;

    PetscDM dmSoln = solution.dmMesh();
    PetscDM dmAux = _auxField->dmMesh();
    PetscDMLabel dmLabel;

    // Pointwise function have been set in DS
    err = DMGetDS(dmSoln, &prob); PYLITH_CHECK_ERROR(err);

    // Get auxiliary data
    err = PetscObjectCompose((PetscObject) dmSoln, "dmAux", (PetscObject) dmAux); PYLITH_CHECK_ERROR(err);
    err = PetscObjectCompose((PetscObject) dmSoln, "A", (PetscObject) _auxField->localVector()); PYLITH_CHECK_ERROR(err);

    // Compute the local residual
    assert(solution.localVector());
    assert(residual->localVector());
    err = DMGetLabel(dmSoln, _label.c_str(), &dmLabel); PYLITH_CHECK_ERROR(err);
    const int labelId = 1;
    const topology::Field::SubfieldInfo& info = solution.subfieldInfo(_field.c_str());

    solution.mesh().view(":mesh.txt:ascii_info_detail"); // :DEBUG:

    PYLITH_COMPONENT_DEBUG("DMPlexComputeBdResidualSingle() for boundary '"<<label()<<"')");
    err = DMPlexComputeBdResidualSingle(dmSoln, t, dmLabel, 1, &labelId, info.index, solution.localVector(), solutionDot.localVector(), residual->localVector()); PYLITH_CHECK_ERROR(err);

    PYLITH_METHOD_END;
} // computeRHSResidual

// ----------------------------------------------------------------------
// Compute RHS Jacobian and preconditioner for G(t,s).
void
pylith::bc::Neumann::computeRHSJacobian(PetscMat jacobianMat,
                                        PetscMat preconMat,
                                        const PylithReal t,
                                        const PylithReal dt,
                                        const pylith::topology::Field& solution) {}

// ----------------------------------------------------------------------
// Compute LHS residual for F(t,s,\dot{s}).
void
pylith::bc::Neumann::computeLHSResidual(pylith::topology::Field* residual,
                                        const PylithReal t,
                                        const PylithReal dt,
                                        const pylith::topology::Field& solution,
                                        const pylith::topology::Field& solutionDot) {}

// ----------------------------------------------------------------------
// Compute LHS Jacobian and preconditioner for F(t,s,\dot{s}) with implicit time-stepping.
void
pylith::bc::Neumann::computeLHSJacobianImplicit(PetscMat jacobianMat,
                                                PetscMat precondMat,
                                                const PylithReal t,
                                                const PylithReal dt,
                                                const PylithReal tshift,
                                                const pylith::topology::Field& solution,
                                                const pylith::topology::Field& solutionDot) {}


// ----------------------------------------------------------------------
// Compute inverse of lumped LHS Jacobian for F(t,s,\dot{s}) with explicit time-stepping.
void
pylith::bc::Neumann::computeLHSJacobianLumpedInv(pylith::topology::Field* jacobianInv,
                                                 const PylithReal t,
                                                 const PylithReal dt,
                                                 const PylithReal tshift,
                                                 const pylith::topology::Field& solution) {}


// ----------------------------------------------------------------------
// Set constants used in finite-element integrations.
void
pylith::bc::Neumann::_setFEConstants(const pylith::topology::Field& solution,
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


// End of file
