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

#include "AbsorbingDampers.hh" // implementation of object methods

#include "AbsorbingDampersAuxiliaryFactory.hh" // USES AuxiliaryFactory

#include "pylith/fekernels/AbsorbingDampers.hh" // USES AbsorbingDampers kernels
#include "pylith/topology/Mesh.hh" // USES Mesh
#include "pylith/topology/Field.hh" // USES Field
#include "pylith/topology/VisitorMesh.hh" // USES VecVisitorMesh
#include "pylith/topology/FieldQuery.hh" // HOLDSA FieldQuery
#include "pylith/topology/CoordsVisitor.hh" // USES CoordsVisitor
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
pylith::bc::AbsorbingDampers::AbsorbingDampers(void) :
    _boundaryMesh(NULL),
    _auxAbsorbingDampersFactory(new pylith::bc::AbsorbingDampersAuxiliaryFactory)
{ // constructor
    _field = "velocity";
} // constructor

// ----------------------------------------------------------------------
// Destructor.
pylith::bc::AbsorbingDampers::~AbsorbingDampers(void) {
    deallocate();
} // destructor

// ----------------------------------------------------------------------
// Deallocate PETSc and local data structures.
void
pylith::bc::AbsorbingDampers::deallocate(void) {
    PYLITH_METHOD_BEGIN;

    IntegratorPointwise::deallocate();

    delete _boundaryMesh; _boundaryMesh = NULL;
    delete _auxAbsorbingDampersFactory; _auxAbsorbingDampersFactory = NULL;

    PYLITH_METHOD_END;
} // deallocate

// ----------------------------------------------------------------------
// Verify configuration is acceptable.
void
pylith::bc::AbsorbingDampers::verifyConfiguration(const pylith::topology::Field& solution) const {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("verifyConfiguration(solution="<<solution.label()<<")");

    if (!solution.hasSubfield(_field.c_str())) {
        std::ostringstream msg;
        msg << "Cannot apply absorbing boundary condition to field '"<< _field
            << "'; field is not in solution.";
        throw std::runtime_error(msg.str());
    } // if

    const pylith::topology::Field::SubfieldInfo& info = solution.subfieldInfo(_field.c_str());
    if (pylith::topology::Field::VECTOR != info.description.vectorFieldType) {
        std::ostringstream msg;
        msg << "Absorbing boundary condition cannot be applied to non-vector field '"<< _field
            << "' in solution.";
        throw std::runtime_error(msg.str());
    } // if

    PYLITH_METHOD_END;
} // verifyConfiguration


// ----------------------------------------------------------------------
// Initialize boundary condition.
void
pylith::bc::AbsorbingDampers::initialize(const pylith::topology::Field& solution) {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("initialize(solution="<<solution.label()<<")");

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

    //_auxField->view("AUXILIARY FIELD"); // :DEBUG:

    PYLITH_METHOD_END;
} // initialize

// ----------------------------------------------------------------------
// Compute RHS residual for G(t,s).
void
pylith::bc::AbsorbingDampers::computeRHSResidual(pylith::topology::Field* residual,
                                                 const PylithReal t,
                                                 const PylithReal dt,
                                                 const pylith::topology::Field& solution) {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("computeRHSResidual(residual="<<residual<<", t="<<t<<", dt="<<dt<<", solution="<<solution.label()<<")");

    _setFEKernelsRHSResidual(solution);

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

    //solution.mesh().view(":mesh.txt:ascii_info_detail"); // :DEBUG:

    PYLITH_COMPONENT_DEBUG("DMPlexComputeBdResidualSingle() for boundary '"<<label()<<"')");
    err = DMPlexComputeBdResidualSingle(dmSoln, t, dmLabel, 1, &labelId, info.index, solution.localVector(), solutionDot.localVector(), residual->localVector()); PYLITH_CHECK_ERROR(err);

    PYLITH_METHOD_END;
} // computeRHSResidual

// ----------------------------------------------------------------------
// Compute RHS Jacobian and preconditioner for G(t,s).
void
pylith::bc::AbsorbingDampers::computeRHSJacobian(PetscMat jacobianMat,
                                                 PetscMat preconMat,
                                                 const PylithReal t,
                                                 const PylithReal dt,
                                                 const pylith::topology::Field& solution) {}

// ----------------------------------------------------------------------
// Compute LHS residual for F(t,s,\dot{s}).
void
pylith::bc::AbsorbingDampers::computeLHSResidual(pylith::topology::Field* residual,
                                                 const PylithReal t,
                                                 const PylithReal dt,
                                                 const pylith::topology::Field& solution,
                                                 const pylith::topology::Field& solutionDot) {}

// ----------------------------------------------------------------------
// Compute LHS Jacobian and preconditioner for F(t,s,\dot{s}) with implicit time-stepping.
void
pylith::bc::AbsorbingDampers::computeLHSJacobianImplicit(PetscMat jacobianMat,
                                                         PetscMat precondMat,
                                                         const PylithReal t,
                                                         const PylithReal dt,
                                                         const PylithReal tshift,
                                                         const pylith::topology::Field& solution,
                                                         const pylith::topology::Field& solutionDot) {}


// ----------------------------------------------------------------------
// Compute inverse of lumped LHS Jacobian for F(t,s,\dot{s}) with explicit time-stepping.
void
pylith::bc::AbsorbingDampers::computeLHSJacobianLumpedInv(pylith::topology::Field* jacobianInv,
                                                          const PylithReal t,
                                                          const PylithReal dt,
                                                          const pylith::topology::Field& solution) {}


// ----------------------------------------------------------------------
// Setup auxiliary subfields (discretization and query fns).
void
pylith::bc::AbsorbingDampers::_auxFieldSetup(const pylith::topology::Field& solution)
{ // _auxFieldSetup
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("_auxFieldSetup(solution="<<solution.label()<<")");

    assert(_auxAbsorbingDampersFactory);
    assert(_normalizer);
    _auxAbsorbingDampersFactory->initialize(_auxField, *_normalizer, solution.spaceDim(),
                                            &solution.subfieldInfo(_field.c_str()).description);

    // :ATTENTION: The order of the factory methods must match the order of the auxiliary subfields in the FE kernels.

    _auxAbsorbingDampersFactory->density();
    _auxAbsorbingDampersFactory->vp();
    _auxAbsorbingDampersFactory->vs();

    PYLITH_METHOD_END;
} // _auxFieldSetup


// ----------------------------------------------------------------------
// Set kernels for RHS residual G(t,s).
void
pylith::bc::AbsorbingDampers::_setFEKernelsRHSResidual(const pylith::topology::Field& solution) const
{ // _setFEKernelsResidual
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("_setFEKernelsResidual(solution="<<solution.label()<<")");

    const PetscDM dmSoln = solution.dmMesh(); assert(dmSoln);
    PetscDS prob = NULL;
    PetscErrorCode err = DMGetDS(dmSoln, &prob); PYLITH_CHECK_ERROR(err);

    PetscBdPointFunc g0 = pylith::fekernels::AbsorbingDampers::g0;
    PetscBdPointFunc g1 = NULL;

    const int fieldIndex = solution.subfieldInfo(_field.c_str()).index;
    err = PetscDSSetBdResidual(prob, fieldIndex, g0, g1); PYLITH_CHECK_ERROR(err);

    PYLITH_METHOD_END;
} // _setFEKernelsResidual


// ----------------------------------------------------------------------
// Get factory for setting up auxliary fields.
pylith::feassemble::AuxiliaryFactory*
pylith::bc::AbsorbingDampers::_auxFactory(void) {
    return _auxAbsorbingDampersFactory;
} // _auxFactory


// End of file
