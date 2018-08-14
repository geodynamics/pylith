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

#include "FaultCohesiveStub.hh" \
    // implementation of object methods

// ----------------------------------------------------------------------
// Default constructor.
pylith::faults::FaultCohesiveStub::FaultCohesiveStub(void)
{}

// ----------------------------------------------------------------------
// Destructor.
pylith::faults::FaultCohesiveStub::~FaultCohesiveStub(void) {
    deallocate();
} // destructor

// ----------------------------------------------------------------------
// Compute RHS residual for G(t,s).
void
pylith::faults::FaultCohesiveStub::computeRHSResidual(pylith::topology::Field* residual,
                                                      const PylithReal t,
                                                      const PylithReal dt,
                                                      const pylith::topology::Field& solution) {}

// ----------------------------------------------------------------------
// Compute RHS Jacobian and preconditioner for G(t,s).
void
pylith::faults::FaultCohesiveStub::computeRHSJacobian(PetscMat jacobianMat,
                                                      PetscMat preconMat,
                                                      const PylithReal t,
                                                      const PylithReal dt,
                                                      const pylith::topology::Field& solution) {}

// ----------------------------------------------------------------------
// Compute LHS residual for F(t,s,\dot{s}).
void
pylith::faults::FaultCohesiveStub::computeLHSResidual(pylith::topology::Field* residual,
                                                      const PylithReal t,
                                                      const PylithReal dt,
                                                      const pylith::topology::Field& solution,
                                                      const pylith::topology::Field& solutionDot) {}

// ----------------------------------------------------------------------
// Compute LHS Jacobian and preconditioner for F(t,s,\dot{s}) with implicit time-stepping.
void
pylith::faults::FaultCohesiveStub::computeLHSJacobianImplicit(PetscMat jacobianMat,
                                                              PetscMat precondMat,
                                                              const PylithReal t,
                                                              const PylithReal dt,
                                                              const PylithReal s_tshift,
                                                              const pylith::topology::Field& solution,
                                                              const pylith::topology::Field& solutionDot) {}


// ----------------------------------------------------------------------
// Compute inverse of lumped LHS Jacobian for F(t,s,\dot{s}) with explicit time-stepping.
void
pylith::faults::FaultCohesiveStub::computeLHSJacobianLumpedInv(pylith::topology::Field* jacobianInv,
                                                               const PylithReal t,
                                                               const PylithReal dt,
                                                               const PylithReal s_tshift,
                                                               const pylith::topology::Field& solution) {}

// ----------------------------------------------------------------------
// Get factory for setting up auxliary fields.
pylith::feassemble::AuxiliaryFactory*
pylith::faults::FaultCohesiveStub::_auxFactory(void) {
    return NULL;
} // _auxFactory


// ----------------------------------------------------------------------
// Setup auxiliary subfields (discretization and query fns).
void
pylith::faults::FaultCohesiveStub::_auxFieldSetup(void) {}


// End of file
