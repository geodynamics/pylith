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
// Copyright (c) 2010-2017 University of California, Davis
//
// See COPYING for license information.
//
// ----------------------------------------------------------------------
//

#include <portinfo>

#include "FaultCohesiveKin.hh" // implementation of object methods

#include "pylith/faults/KinSrc.hh" // USES KinSrc
#include "pylith/faults/AuxiliaryFactory.hh" // USES AuxiliaryFactory

#include "pylith/topology/Mesh.hh" // USES Mesh
#include "pylith/topology/Field.hh" // USES Field

#include "pylith/fekernels/FaultCohesiveKin.hh" // USES FaultCohesiveKin

#include "pylith/utils/EventLogger.hh" // USES EventLogger
#include "pylith/utils/journals.hh" // USES PYLITH_COMPONENT_*

#include "spatialdata/geocoords/CoordSys.hh" // USES CoordSys
#include "spatialdata/spatialdb/SpatialDB.hh" // USES SpatialDB



#include <cmath> // USES pow(), sqrt()
#include <strings.h> // USES strcasecmp()
#include <cstring> // USES strlen()
#include <cstdlib> // USES atoi()
#include <cassert> // USES assert()
#include <sstream> // USES std::ostringstream
#include <stdexcept> \
    // USES std::runtime_error

//#define DETAILED_EVENT_LOGGING

// ----------------------------------------------------------------------
// Default constructor.
pylith::faults::FaultCohesiveKin::FaultCohesiveKin(void)
{ // constructor
} // constructor

// ----------------------------------------------------------------------
// Destructor.
pylith::faults::FaultCohesiveKin::~FaultCohesiveKin(void)
{ // destructor
    deallocate();
} // destructor

// ----------------------------------------------------------------------
// Deallocate PETSc and local data structures.
void
pylith::faults::FaultCohesiveKin::deallocate(void) {
    FaultCohesive::deallocate();

    _eqSrcs.clear(); // :TODO: Use shared pointers for earthquake sources
} // deallocate

// ----------------------------------------------------------------------
// Set kinematic earthquake source.
void
pylith::faults::FaultCohesiveKin::eqsrcs(const char* const * names,
                                         const int numNames,
                                         KinSrc** sources,
                                         const int numSources) {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("eqsrcs(names="<<names<<", numNames="<<numNames<<", sources="<<sources<<", numSources="<<numSources<<")");

    assert(numNames == numSources);

    // :TODO: Use shared pointers for earthquake sources
    _eqSrcs.clear();
    for (int i = 0; i < numSources; ++i) {
        if (0 == sources[i]) {
            throw std::runtime_error("Null earthquake source.");
        }
        _eqSrcs[std::string(names[i])] = sources[i];
    } // for

    PYLITH_METHOD_END;
} // eqsrcs


// ----------------------------------------------------------------------
// Initialize fault. Determine orientation and setup boundary
void
pylith::faults::FaultCohesiveKin::initialize(const pylith::topology::Field& solution) {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("initialize(solution="<<solution.label()<<")");

    assert(_normalizer);

    FaultCohesive::initialize(solution);

    assert(_auxField);
    const srcs_type::const_iterator srcsEnd = _eqSrcs.end();
    for (srcs_type::iterator s_iter = _eqSrcs.begin(); s_iter != srcsEnd; ++s_iter) {
        KinSrc* src = s_iter->second;
        assert(src);
        src->initialize(*_auxField, *_normalizer, solution.mesh().coordsys());
    } // for

    PYLITH_METHOD_END;
} // initialize

// ----------------------------------------------------------------------
// Update auxiliary fields at beginning of time step.
void
pylith::faults::FaultCohesiveKin::prestep(const double t,
                                          const double dt) {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("prestep(t="<<t<<", dt="<<dt<<")");

    // Update slip subfield.

    // Compute slip field at current time step
    const srcs_type::const_iterator srcsEnd = _eqSrcs.end();
    for (srcs_type::iterator s_iter = _eqSrcs.begin(); s_iter != srcsEnd; ++s_iter) {
        KinSrc* src = s_iter->second; assert(src);
        src->slip(_auxField, t);
    } // for

    PYLITH_METHOD_END;
} // prestep


// ----------------------------------------------------------------------
// Compute RHS residual for G(t,s).
void
pylith::faults::FaultCohesiveKin::computeRHSResidual(pylith::topology::Field* residual,
                                                     const PylithReal t,
                                                     const PylithReal dt,
                                                     const pylith::topology::Field& solution) {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("computeRHSResidual(residual="<<residual<<", t="<<t<<", dt="<<dt<<", solution="<<solution.label()<<")");

    assert(residual);
    assert(_auxField);

    pylith::topology::Field solutionDot(solution.mesh()); // No dependence on time derivative of solution in RHS.
    solutionDot.label("solution_dot");

    const int i_lagrange = solution.subfieldInfo("lagrange_multiplier_fault").index;
    const int i_dispvel = (solution.hasSubfield("velocity")) ?
                          solution.subfieldInfo("velocity").index : solution.subfieldInfo("displacement").index;

    PetscErrorCode err;
    PetscDM dmSoln = solution.dmMesh();
    PetscDM dmAux = _auxField->dmMesh();

    // Pointwise functions have been set in DS
    PetscDS prob = NULL;
    err = DMGetDS(dmSoln, &prob); PYLITH_CHECK_ERROR(err);

    // Get auxiliary data
    err = PetscObjectCompose((PetscObject) dmSoln, "dmAux", (PetscObject) dmAux); PYLITH_CHECK_ERROR(err);
    err = PetscObjectCompose((PetscObject) dmSoln, "A", (PetscObject) _auxField->localVector()); PYLITH_CHECK_ERROR(err);

    _setFEConstants(solution, dt);

    // Compute the local residual
    assert(solution.localVector());
    assert(residual->localVector());

    PetscDMLabel dmLabel = NULL;
    const int labelId = 1;
    std::string faultLabel;

    PYLITH_COMPONENT_DEBUG("DMPlexComputeBdResidualSingle() on the positive side of fault '"<<label()<<"'.");
    _setFEKernelsRHSResidualFaultPositive(solution);
    faultLabel = std::string(label()) + std::string("_fault_positive");
    err = DMGetLabel(dmSoln, faultLabel.c_str(), &dmLabel); PYLITH_CHECK_ERROR(err);
    err = DMPlexComputeBdResidualSingle(dmSoln, t, dmLabel, 1, &labelId, i_dispvel, solution.localVector(),
                                        solutionDot.localVector(), residual->localVector()); PYLITH_CHECK_ERROR(err);
    err = DMPlexComputeBdResidualSingle(dmSoln, t, dmLabel, 1, &labelId, i_lagrange, solution.localVector(),
                                        solutionDot.localVector(), residual->localVector()); PYLITH_CHECK_ERROR(err);

    PYLITH_COMPONENT_DEBUG("DMPlexComputeBdResidualSingle() on the negative side of fault '"<<label()<<"'.");
    _setFEKernelsRHSResidualFaultNegative(solution);
    faultLabel = std::string(label()) + std::string("_fault_negative");
    err = DMGetLabel(dmSoln, faultLabel.c_str(), &dmLabel); PYLITH_CHECK_ERROR(err);
    err = DMPlexComputeBdResidualSingle(dmSoln, t, dmLabel, 1, &labelId, i_dispvel, solution.localVector(),
                                        solutionDot.localVector(), residual->localVector()); PYLITH_CHECK_ERROR(err);
    err = DMPlexComputeBdResidualSingle(dmSoln, t, dmLabel, 1, &labelId, i_lagrange, solution.localVector(),
                                        solutionDot.localVector(), residual->localVector()); PYLITH_CHECK_ERROR(err);

    PYLITH_METHOD_END;
} // computeRHSResidual


// ----------------------------------------------------------------------
// Compute RHS Jacobian and preconditioner for G(t,s).
void
pylith::faults::FaultCohesiveKin::computeRHSJacobian(PetscMat jacobianMat,
                                                     PetscMat precondMat,
                                                     const PylithReal t,
                                                     const PylithReal dt,
                                                     const pylith::topology::Field& solution) {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("computeRHSJacobian(jacobianMat="<<jacobianMat<<", precondMat="<<precondMat<<", t="<<t<<", dt="<<dt<<", solution="<<solution.label()<<")");

    assert(jacobianMat);
    assert(precondMat);
    assert(_auxField);

    pylith::topology::Field solutionDot(solution.mesh()); // No dependence on time derivative of solution in RHS.
    solutionDot.label("solution_dot");
    const PylithReal tshift = 0.0; // No dependence on time derivative of solution in RHS, so shift isn't applicable.

    PetscErrorCode err;
    PetscDM dmSoln = solution.dmMesh();
    PetscDM dmAux = _auxField->dmMesh();

    // Pointwise function have been set in DS
    PetscDS prob = NULL;
    err = DMGetDS(dmSoln, &prob); PYLITH_CHECK_ERROR(err);

    // Get auxiliary data
    err = PetscObjectCompose((PetscObject) dmSoln, "dmAux", (PetscObject) dmAux); PYLITH_CHECK_ERROR(err);
    err = PetscObjectCompose((PetscObject) dmSoln, "A", (PetscObject) auxField()->localVector()); PYLITH_CHECK_ERROR(err);

    _setFEConstants(solution, dt);

    // Compute the local Jacobian
    assert(solution.localVector());

    PetscDMLabel dmLabel = NULL;
    const int labelId = 1;
    std::string faultLabel;

#if 0
    PYLITH_COMPONENT_DEBUG("DMPlexComputeBdJacobianSingle() on the positive side of fault ''"<<label()<<"')");
    _setFEKernelsRHSJacobianFaultPositive(solution);
    faultLabel = std::string(label()) + std::string("_fault_positive");
    err = DMGetLabel(dmSoln, faultLabel.c_str(), &dmLabel); PYLITH_CHECK_ERROR(err);
    err = DMPlexComputeBdJacobianSingle(dmSoln, t, tshift, dmLabel, 1, &labelId, solution.localVector(),
                                        solutionDot.localVector(), jacobianMat, precondMat, NULL); PYLITH_CHECK_ERROR(err);

    PYLITH_COMPONENT_DEBUG("DMPlexComputeBdJacobianSingle() on the negative side of fault ''"<<label()<<"')");
    _setFEKernelsRHSJacobianFaultNegative(solution);
    faultLabel = std::string(label()) + std::string("_fault_negative");
    err = DMGetLabel(dmSoln, faultLabel.c_str(), &dmLabel); PYLITH_CHECK_ERROR(err);
    err = DMPlexComputeBdJacobianSingle(dmSoln, t, tshift, dmLabel, 1, &labelId, solution.localVector(),
                                        solutionDot.localVector(), jacobianMat, precondMat, NULL); PYLITH_CHECK_ERROR(err);

#else
    PYLITH_COMPONENT_ERROR(":TODO: @matt @brad Waiting for Matt to implement DMPlexComputeBdJacobianSingle().");
#endif

    _needNewRHSJacobian = false;

    PYLITH_METHOD_END;
} // computeRHSJacobian


// ----------------------------------------------------------------------
// Compute LHS residual for F(t,s,\dot{s}).
void
pylith::faults::FaultCohesiveKin::computeLHSResidual(pylith::topology::Field* residual,
                                                     const PylithReal t,
                                                     const PylithReal dt,
                                                     const pylith::topology::Field& solution,
                                                     const pylith::topology::Field& solutionDot) {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("computeLHSResidual(residual="<<residual<<", t="<<t<<", dt="<<dt<<", solution="<<solution.label()<<", solutionDot="<<solutionDot.label()<<") empty method");

    // No contribution to LHS residual.

    PYLITH_METHOD_END;
} // computeLHSResidual


// ----------------------------------------------------------------------
// Compute LHS Jacobian and preconditioner for F(t,s,\dot{s}) with implicit time-stepping.
void
pylith::faults::FaultCohesiveKin::computeLHSJacobianImplicit(PetscMat jacobianMat,
                                                             PetscMat precondMat,
                                                             const PylithReal t,
                                                             const PylithReal dt,
                                                             const PylithReal tshift,
                                                             const pylith::topology::Field& solution,
                                                             const pylith::topology::Field& solutionDot) {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG(
        "computeLHSJacobianImplicit(jacobianMat="<<jacobianMat<<", precondMat="<<precondMat<<", t="<<t<<", dt="<<dt<<", tshift="<<tshift<<", solution="<<solution.label()<<", solutionDot="<<solutionDot.label()<<") empty method");

    // No contribution to LHS Jacobian.

    PYLITH_METHOD_END;
} // computeLHSJacobianImplicit


// ----------------------------------------------------------------------
// Compute inverse of lumped LHS Jacobian for F(t,s,\dot{s}) with explicit time-stepping.
void
pylith::faults::FaultCohesiveKin::computeLHSJacobianLumpedInv(pylith::topology::Field* jacobianInv,
                                                              const PylithReal t,
                                                              const PylithReal dt,
                                                              const PylithReal tshift,
                                                              const pylith::topology::Field& solution) {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("computeLHSJacobianLumpedInv(jacobianInv="<<jacobianInv<<", t="<<t<<", dt="<<dt<<", solution="<<solution.label()<<") empty method");

    // No contribution to LHS Jacobian.

    PYLITH_METHOD_END;
} // computeLHSJacobianLumpedInv


// ----------------------------------------------------------------------
// Setup auxiliary subfields (discretization and query fns).
void
pylith::faults::FaultCohesiveKin::_auxFieldSetup(void) {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("_auxFieldSetup()");

    assert(_auxFaultFactory);
    assert(_normalizer);
    assert(_faultMesh);

    const spatialdata::geocoords::CoordSys* cs = _faultMesh->coordsys();assert(cs);
    const int spaceDim = cs->spaceDim();
    _auxFaultFactory->initialize(_auxField, *_normalizer, spaceDim);

    // :ATTENTION: The order for adding subfields must match the order of the auxiliary fields in the FE kernels.

    _auxFaultFactory->normalDir(); // 0
    _auxFaultFactory->strikeDir(); // 1
    if (3 == spaceDim) {
        _auxFaultFactory->upDipDir(); // 2
    } // if
    _auxFaultFactory->slip(); // numA-1

    PYLITH_METHOD_END;
} // _auxFieldSetup


// ----------------------------------------------------------------------
// Set pointwise functions for computing contributions of the positive
// side of the fault to the RHS residual.
void
pylith::faults::FaultCohesiveKin::_setFEKernelsRHSResidualFaultPositive(const pylith::topology::Field& solution) {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("_setFEKernelsRHSResidualFaultPositive(solution="<<solution.label()<<")");

    const PetscDM dm = solution.dmMesh(); assert(dm);
    PetscDS prob = NULL;
    PetscErrorCode err = DMGetDS(dm, &prob); PYLITH_CHECK_ERROR(err);

    // Elasticity equation (displacement/velocity).
    const PetscInt i_dispvel = (solution.hasSubfield("velocity")) ?
                               solution.subfieldInfo("velocity").index : solution.subfieldInfo("displacement").index;
    const PetscPointFunc g0u = pylith::fekernels::FaultCohesiveKin::g0u_pos;
    const PetscPointFunc g1u = NULL;
    err = PetscDSSetResidual(prob, i_dispvel, g0u, g1u); PYLITH_CHECK_ERROR(err);

    // Fault slip constraint equation.
    const PetscInt i_lagrange = solution.subfieldInfo("lagrange_multiplier_fault").index;
    const PetscPointFunc g0l = pylith::fekernels::FaultCohesiveKin::g0l_pos;
    const PetscPointFunc g1l = NULL;
    err = PetscDSSetResidual(prob, i_lagrange,  g0l, g1l); PYLITH_CHECK_ERROR(err);

    PYLITH_METHOD_END;
} // _setFEKernelsRHSResidualFaultPositive

// ----------------------------------------------------------------------
// Set pointwise functions for computing contributions of the positive
// side of the fault to the RHS residual.
void
pylith::faults::FaultCohesiveKin::_setFEKernelsRHSResidualFaultNegative(const pylith::topology::Field& solution) {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("_setFEKernelsRHSResidualFaultNegative(solution="<<solution.label()<<")");

    const PetscDM dm = solution.dmMesh(); assert(dm);
    PetscDS prob = NULL;
    PetscErrorCode err = DMGetDS(dm, &prob); PYLITH_CHECK_ERROR(err);

    // Elasticity equation (displacement/velocity).
    const PetscInt i_dispvel = (solution.hasSubfield("velocity")) ?
                               solution.subfieldInfo("velocity").index : solution.subfieldInfo("displacement").index;
    const PetscPointFunc g0u = pylith::fekernels::FaultCohesiveKin::g0u_neg;
    const PetscPointFunc g1u = NULL;
    err = PetscDSSetResidual(prob, i_dispvel, g0u, g1u); PYLITH_CHECK_ERROR(err);

    // Fault slip constraint equation.
    const PetscInt i_lagrange = solution.subfieldInfo("lagrange_multiplier_fault").index;
    const PetscPointFunc g0l = pylith::fekernels::FaultCohesiveKin::g0l_neg;
    const PetscPointFunc g1l = NULL;
    err = PetscDSSetResidual(prob, i_lagrange,  g0l, g1l); PYLITH_CHECK_ERROR(err);

    PYLITH_METHOD_END;
} // _setFEKernelsRHSResidualFaultPositive

// ----------------------------------------------------------------------
// Set pointwise functions for computing contributions of the positive
// side of the fault to the RHS Jacobian.
void
pylith::faults::FaultCohesiveKin::_setFEKernelsRHSJacobianFaultPositive(const pylith::topology::Field& solution) {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("_setFEKernelsRHSJacobianFaultPositive(solution="<<solution.label()<<")");

    const PetscDM dm = solution.dmMesh(); assert(dm);
    PetscDS prob = NULL;
    PetscErrorCode err = DMGetDS(dm, &prob); PYLITH_CHECK_ERROR(err);

    const PetscInt i_dispvel = (solution.hasSubfield("velocity")) ?
                               solution.subfieldInfo("velocity").index : solution.subfieldInfo("displacement").index;
    const PetscInt i_lagrange = solution.subfieldInfo("lagrange_multiplier_fault").index;

    const PetscPointJac Jg0ul = pylith::fekernels::FaultCohesiveKin::Jg0ul_pos;
    const PetscPointJac Jg1ul = NULL;
    const PetscPointJac Jg2ul = NULL;
    const PetscPointJac Jg3ul = NULL;

    const PetscPointJac Jg0lu = pylith::fekernels::FaultCohesiveKin::Jg0lu_pos;
    const PetscPointJac Jg1lu = NULL;
    const PetscPointJac Jg2lu = NULL;
    const PetscPointJac Jg3lu = NULL;

    err = PetscDSSetJacobian(prob, i_dispvel, i_lagrange, Jg0ul, Jg1ul, Jg2ul, Jg3ul); PYLITH_CHECK_ERROR(err);
    err = PetscDSSetJacobian(prob, i_lagrange, i_dispvel, Jg0lu, Jg1lu, Jg2lu, Jg3lu); PYLITH_CHECK_ERROR(err);

    PYLITH_METHOD_END;
} // _setFEKernelsRHSJacobianFaultPositive


// ----------------------------------------------------------------------
// Set pointwise functions for computing contributions of the positive
// side of the fault to the RHS Jacobian.
void
pylith::faults::FaultCohesiveKin::_setFEKernelsRHSJacobianFaultNegative(const pylith::topology::Field& solution) {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("_setFEKernelsRHSJacobianFaultNegative(solution="<<solution.label()<<")");

    const PetscDM dm = solution.dmMesh(); assert(dm);
    PetscDS prob = NULL;
    PetscErrorCode err = DMGetDS(dm, &prob); PYLITH_CHECK_ERROR(err);

    const PetscInt i_dispvel = (solution.hasSubfield("velocity")) ?
                               solution.subfieldInfo("velocity").index : solution.subfieldInfo("displacement").index;
    const PetscInt i_lagrange = solution.subfieldInfo("lagrange_multiplier_fault").index;

    const PetscPointJac Jg0ul = pylith::fekernels::FaultCohesiveKin::Jg0ul_neg;
    const PetscPointJac Jg1ul = NULL;
    const PetscPointJac Jg2ul = NULL;
    const PetscPointJac Jg3ul = NULL;

    const PetscPointJac Jg0lu = pylith::fekernels::FaultCohesiveKin::Jg0lu_neg;
    const PetscPointJac Jg1lu = NULL;
    const PetscPointJac Jg2lu = NULL;
    const PetscPointJac Jg3lu = NULL;

    err = PetscDSSetJacobian(prob, i_dispvel, i_lagrange, Jg0ul, Jg1ul, Jg2ul, Jg3ul); PYLITH_CHECK_ERROR(err);
    err = PetscDSSetJacobian(prob, i_lagrange, i_dispvel, Jg0lu, Jg1lu, Jg2lu, Jg3lu); PYLITH_CHECK_ERROR(err);

    PYLITH_METHOD_END;
} // _setFEKernelsRHSJacobianFaultPositive

// End of file
