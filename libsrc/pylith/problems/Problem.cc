// -*- C++ -*-
//
// ======================================================================
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
// ======================================================================
//

#include <portinfo>

#include "Problem.hh" // implementation of class methods

#include "pylith/topology/Mesh.hh" // USES Mesh
#include "pylith/topology/Field.hh" // USES Field

#include "pylith/feassemble/IntegratorPointwise.hh" // USES IntegratorPointwise
#include "pylith/topology/Jacobian.hh" // USES Jacobian

#include "pylith/utils/error.hh" // USES PYLITH_CHECK_ERROR

#include "journal/debug.h" // USES journal::debug_t

#include <cassert> // USES assert()

// ----------------------------------------------------------------------
// Constructor
pylith::problems::Problem::Problem(void) :
  _solution(0),
  _jacobianLHSLumpedInv(0),
  _integrators(0),
  _constraints(0)
{ // constructor
} // constructor

// ----------------------------------------------------------------------
// Destructor
pylith::problems::Problem::~Problem(void)
{ // destructor
  deallocate();
} // destructor

// ----------------------------------------------------------------------
// Deallocate PETSc and local data structures.
void
pylith::problems::Problem::deallocate(void)
{ // deallocate
  PYLITH_METHOD_BEGIN;

  _solution = 0; // Held by Python. :KLUDGE: :TODO: Use shared pointer.
  delete _jacobianLHSLumpedInv; _jacobianLHSLumpedInv = 0;

  PYLITH_METHOD_END;
} // deallocate
  
// ----------------------------------------------------------------------
// Set integrators over the mesh.
void
pylith::problems::Problem::integrators(pylith::feassemble::IntegratorPointwise* integratorArray[],
				       const int numIntegrators)
{ // integrators
  assert( (!integratorArray && 0 == numIntegrators) || (integratorArray && 0 < numIntegrators) );

  _integrators.resize(numIntegrators);
  for (int i=0; i < numIntegrators; ++i) {
    _integrators[i] = integratorArray[i];
  } // for
} // integrators
  
// ----------------------------------------------------------------------
// Set constraints over the mesh.
void
pylith::problems::Problem::constraints(pylith::feassemble::Constraint* constraintArray[],
				       const int numConstraints)
{ // constraints
  assert( (!constraintArray && 0 == numConstraints) || (constraintArray && 0 < numConstraints) );

  _constraints.resize(numConstraints);
  for (int i=0; i < numConstraints; ++i) {
    _constraints[i] = constraintArray[i];
  } // for
} // constraints
  
// ----------------------------------------------------------------------
// Initialize.
void
pylith::problems::Problem::initialize(void)
{ // initialize
} // initialize

// ----------------------------------------------------------------------
// Compute RHS residual for G(t,s).
void
pylith::problems::Problem::computeRHSResidual(PetscVec residualVec,
					      const PylithReal t,
					      const PylithReal dt,
					      PetscVec solutionVec)
{ // computeRHSResidual
  PYLITH_METHOD_BEGIN;

  journal::debug_t debug("problem");
  debug << journal::at(__HERE__)
	<< "Problem::computeRHSResidual(t="<<t<<", dt="<<dt<<", solutionVec="<<solutionVec<<", residualVec="<<residualVec<<")" << journal::endl;

  assert(residualVec);
  assert(solutionVec);
  assert(_solution);

  // Update PyLith view of the solution.
  _solution->scatterGlobalToLocal(solutionVec);

  // Sum residual contributions across integrators.
  const size_t numIntegrators = _integrators.size();
  assert(numIntegrators > 0); // must have at least 1 integrator
  for (size_t i=0; i < numIntegrators; ++i) {
    _integrators[i]->computeRHSResidual(residualVec, t, dt, *_solution);
  } // for

  PYLITH_METHOD_END;
} // computeRHSResidual

// ----------------------------------------------------------------------
// Compute RHS Jacobian for G(t,s).
void
pylith::problems::Problem::computeRHSJacobian(PetscMat jacobianMat,
					      PetscMat precondMat,
					      const PylithReal t,
					      const PylithReal dt,
					      PetscVec solutionVec)
{ // computeRHSJacobian
  PYLITH_METHOD_BEGIN;

  journal::debug_t debug("problem");
  debug << journal::at(__HERE__)
	<< "Problem::computeRHSJacobian(t="<<t<<", dt="<<dt<<", solutionVec="<<solutionVec<<", jacobianMat="<<jacobianMat<<", precondMat="<<precondMat<<")" << journal::endl;

  // :KLUDGE: Should add check to see if we need to compute Jacobian

  // Update PyLith view of the solution.
  assert(_solution);
  _solution->scatterGlobalToLocal(solutionVec);

  // Sum Jacobian contributions across integrators.
  const size_t numIntegrators = _integrators.size();
  for (size_t i=0; i < numIntegrators; ++i) {
    _integrators[i]->computeRHSJacobian(jacobianMat, precondMat, t, dt, *_solution);
  } // for
  
  PYLITH_METHOD_END;
} // computeRHSJacobian


// ----------------------------------------------------------------------
// Compute LHS residual for F(t,s,\dot{s}).
void
pylith::problems::Problem::computeLHSResidual(PetscVec residualVec,
					      const PylithReal t,
					      const PylithReal dt,
					      PetscVec solutionVec,
					      PetscVec solutionDotVec)
{ // computeLHSResidual
  PYLITH_METHOD_BEGIN;

  journal::debug_t debug("problem");
  debug << journal::at(__HERE__)
	<< "Problem::computeLHSResidual(t="<<t<<", dt="<<dt<<", solutionVec="<<solutionVec<<", solutionDotVec"<<solutionDotVec<<", residualVec="<<residualVec<<")" << journal::endl;

  assert(residualVec);
  assert(solutionVec);
  assert(solutionDotVec);
  assert(_solution);

  // Update PyLith view of the solution.
  _solution->scatterGlobalToLocal(solutionVec);

  // Sum residual across integrators.
  const int numIntegrators = _integrators.size();
  assert(numIntegrators > 0); // must have at least 1 integrator
  for (int i=0; i < numIntegrators; ++i) {
    _integrators[i]->computeLHSResidual(residualVec, t, dt, *_solution, solutionDotVec);
  } // for

  PYLITH_METHOD_END;
} // computeLHSResidual


// ----------------------------------------------------------------------
// Compute LHS Jacobian for F(t,s,\dot{s}) for implicit time stepping.
void
pylith::problems::Problem::computeLHSJacobianImplicit(PetscMat jacobianMat,
						      PetscMat precondMat,
						      const PylithReal t,
						      const PylithReal dt,
						      const PylithReal tshift,
						      PetscVec solutionVec,
						      PetscVec solutionDotVec)
{ // computeLHSJacobianImplicit
  PYLITH_METHOD_BEGIN;

  journal::debug_t debug("problem");
  debug << journal::at(__HERE__)
	<< "Problem::computeLHSJacobianImplicit(t="<<t<<", dt="<<dt<<", tshift="<<tshift<<", solutionVec="<<solutionVec<<", solutionDotVec="<<solutionDotVec<<", jacobianMat="<<jacobianMat<<", precondMat="<<precondMat<<")" << journal::endl;

  // :KLUDGE: :TODO: Should add check to see if we need to compute Jacobian

  // Update PyLith view of the solution.
  assert(_solution);
  _solution->scatterGlobalToLocal(solutionVec);

  // Sum Jacobian contributions across integrators.
  const size_t numIntegrators = _integrators.size();
  for (size_t i=0; i < numIntegrators; ++i) {
    _integrators[i]->computeLHSJacobianImplicit(jacobianMat, precondMat, t, dt, tshift, *_solution, solutionDotVec);
  } // for
  
  PYLITH_METHOD_END;
} // computeLHSJacobianImplicit

// ----------------------------------------------------------------------
// Compute LHS Jacobian for F(t,s,\dot{s}) for explicit time stepping.
void
pylith::problems::Problem::computeLHSJacobianLumpedInv(const PylithReal t,
						       const PylithReal dt,
						       PetscVec solutionVec)
{ // computeLHSJacobianLumpedInv
  PYLITH_METHOD_BEGIN;

  // :KLUDGE: :TODO: Should add check to see if we need to compute Jacobian

  // Set jacobian to zero.
  assert(_jacobianLHSLumpedInv);
  _jacobianLHSLumpedInv->zeroAll();

  // Update PyLith view of the solution.
  assert(_solution);
  _solution->scatterGlobalToLocal(solutionVec);

  // Sum Jacobian contributions across integrators.
  const size_t numIntegrators = _integrators.size();
  for (size_t i=0; i < numIntegrators; ++i) {
    _integrators[i]->computeLHSJacobianLumpedInv(_jacobianLHSLumpedInv, t, dt, *_solution);
  } // for
  
  PYLITH_METHOD_END;
} // computeLHSJacobianLumpedInv

// End of file 
