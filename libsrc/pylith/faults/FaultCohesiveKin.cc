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

#include "pylith/topology/Mesh.hh" // USES Mesh
#include "pylith/topology/Field.hh" // USES Field

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

    PYLITH_COMPONENT_ERROR(":TODO: @brad Implement computeRHSResidual().");

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

    PYLITH_COMPONENT_ERROR(":TODO: @brad Implement computeRHSResidual().");

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

    // Default is to do nothing.

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

    // Default is to do nothing.

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
    PYLITH_COMPONENT_DEBUG("computeLHSJacobianLumpedInv(jacobianInv="<<jacobianInv<<", t="<<t<<", dt="<<dt<<", solution="<<solution.label()<<")");

    // Default is to do nothing.

    PYLITH_METHOD_END;
} // computeLHSJacobianLumpedInv


// End of file
