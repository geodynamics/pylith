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
#include <cassert> // USES assert()

// ----------------------------------------------------------------------
// Constructor
pylith::problems::Problem::Problem(void) :
  _solution(0),
  _solutionDot(0),
  _residualRHS(0),
  _residualLHS(0),
  _jacobianRHS(0),
  _jacobianLHS(0),
  _preconditionerRHS(0),
  _preconditionerLHS(0),
  _integrators(0),
  _constraints(0),
  _customConstraintPCMat(0)
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
  _solutionDot = 0; // Held by Python. :KLUDGE: :TODO: Use shared pointer.
  delete _residualRHS; _residualRHS = 0;
  delete _residualLHS; _residualLHS = 0;
  delete _jacobianRHS; _jacobianRHS = 0;
  delete _jacobianLHS; _jacobianLHS = 0;
  delete _preconditionerRHS; _preconditionerRHS = 0;
  delete _preconditionerLHS; _preconditionerLHS = 0;

#if 0   // :KLUDGE: Assume Solver deallocates matrix.
  PetscErrorCode err = 0;
  if (_customConstraintPCMat) {
    err = PetscObjectDereference((PetscObject) _customConstraintPCMat);PYLITH_CHECK_ERROR(err);
    _customConstraintPCMat = 0;
  } // if
#else
  _customConstraintPCMat = 0;
#endif

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
// Set handle to preconditioner.
void
pylith::problems::Problem::customPCMatrix(PetscMat& mat)
{ // preconditioner
  _customConstraintPCMat = mat;

#if 0 // :KLUDGE: Assume solver deallocates matrix
  PetscErrorCode err = 0;
  err = PetscObjectReference((PetscObject) mat); PYLITH_CHECK_ERROR(err);
#endif
} // preconditioner

// ----------------------------------------------------------------------
// Initialize.
void
pylith::problems::Problem::initialize(void)
{ // initialize
} // initialize

// ----------------------------------------------------------------------
// Compute RHS residual for G(t,s).
void
pylith::problems::Problem::computeRHSResidual(const PylithReal t,
					      const PylithReal dt,
					      PetscVec solutionVec,
					      PetscVec residualVec)
{ // computeRHSResidual
  PYLITH_METHOD_BEGIN;

  assert(residualVec);
  assert(_residualRHS);
  assert(solutionVec);
  assert(_solution);

  // Update PyLith view of the solution.
  _solution->scatterGlobalToLocal(solutionVec);

  // Set residual to zero.
  _residualRHS->zeroAll();
  
  // Sum residual contributions across integrators.
  const size_t numIntegrators = _integrators.size();
  assert(numIntegrators > 0); // must have at least 1 integrator
  for (size_t i=0; i < numIntegrators; ++i) {
    _integrators[i]->computeRHSResidual(_residualRHS, t, dt, *_solution);
  } // for

  // Update PETSc view of residual
  PetscErrorCode err;
  PetscDM dmMesh = _residualRHS->dmMesh();
  err = DMLocalToGlobalBegin(dmMesh, _residualRHS->localVector(), ADD_VALUES, residualVec);PYLITH_CHECK_ERROR(err);
  err = DMLocalToGlobalEnd(dmMesh, _residualRHS->localVector(), ADD_VALUES, residualVec);PYLITH_CHECK_ERROR(err);

  PYLITH_METHOD_END;
} // computeRHSResidual

// ----------------------------------------------------------------------
// Compute RHS Jacobian for G(t,s).
void
pylith::problems::Problem::computeRHSJacobian(const PylithReal t,
					      const PylithReal dt,
					      PetscVec solutionVec,
					      PetscMat jacobianMat,
					      PetscMat precondMat)
{ // computeRHSJacobian
  PYLITH_METHOD_BEGIN;

  assert(_jacobianRHS);
  assert(_solution);

  // :KLUDGE: Should add check to see if we need to compute Jacobian

  // Update PyLith view of the solution.
  _solution->scatterGlobalToLocal(solutionVec);

  // Set jacobian to zero.
  _jacobianRHS->zero();

  // Sum Jacobian contributions across integrators.
  const size_t numIntegrators = _integrators.size();
  for (size_t i=0; i < numIntegrators; ++i) {
    _integrators[i]->computeRHSJacobian(_jacobianRHS, _preconditionerRHS, t, dt, *_solution);
  } // for
  
  // Assemble jacobian.
  _jacobianRHS->assemble("final_assembly");

#if 0 // FIX THIS
  if (_customConstraintPCMat) {
    // Recalculate preconditioner.
    for (size_t i=0; i < numIntegrators; ++i) {
      _integrators[i]->computeRHSPreconditioner(&_customConstraintPCMat, _jacobianRHS, t, dt, *_solution);
    } // for

    MatAssemblyBegin(_customConstraintPCMat, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(_customConstraintPCMat, MAT_FINAL_ASSEMBLY);

#if 0 // debugging
    std::cout << "Preconditioner Matrix" << std::endl;
    MatView(_customConstraintPCMat, PETSC_VIEWER_STDOUT_WORLD);
#endif
  } // if
#endif

  PYLITH_METHOD_END;
} // computeRHSJacobian


// ----------------------------------------------------------------------
// Compute LHS residual for F(t,s,\dot{s}).
void
pylith::problems::Problem::computeLHSResidual(const PylithReal t,
					      const PylithReal dt,
					      PetscVec solutionVec,
					      PetscVec solutionDotVec,
					      PetscVec residualVec)
{ // computeLHSResidual
  PYLITH_METHOD_BEGIN;

  assert(residualVec);
  assert(_residualLHS);

  assert(solutionVec);
  assert(_solution);
  _solution->scatterGlobalToLocal(solutionVec);

  assert(solutionDotVec);
  assert(_solutionDot);
  _solutionDot->scatterGlobalToLocal(solutionDotVec);

  // Set residual to zero.
  _residualLHS->zeroAll();
  
  // Sum residual across integrators.
  const int numIntegrators = _integrators.size();
  assert(numIntegrators > 0); // must have at least 1 integrator
  for (int i=0; i < numIntegrators; ++i) {
    _integrators[i]->computeLHSResidual(_residualLHS, t, dt, *_solution, *_solutionDot);
  } // for

  // Update PETSc view of residual
  PetscErrorCode err;
  PetscDM dmMesh = _residualRHS->dmMesh();
  err = DMLocalToGlobalBegin(dmMesh, _residualLHS->localVector(), ADD_VALUES, residualVec);PYLITH_CHECK_ERROR(err);
  err = DMLocalToGlobalEnd(dmMesh, _residualLHS->localVector(), ADD_VALUES, residualVec);PYLITH_CHECK_ERROR(err);

  PYLITH_METHOD_END;
} // computeLHSResidual


// ----------------------------------------------------------------------
// Compute LHS Jacobian for F(t,s,\dot{s}) for implicit time stepping.
void
pylith::problems::Problem::computeLHSJacobianImplicit(const PylithReal t,
						      const PylithReal dt,
						      const PylithReal tshift,
						      PetscVec solutionVec,
						      PetscVec solutionDotVec,
						      PetscMat jacobianMat,
						      PetscMat precondMat)
{ // computeLHSJacobianImplicit
  PYLITH_METHOD_BEGIN;

  assert(_jacobianLHS);
  assert(_solution);
  assert(_solutionDot);

  // :KLUDGE: Should add check to see if we need to compute Jacobian

  // Update PyLith view of the solution.
  _solution->scatterGlobalToLocal(solutionVec);
  _solutionDot->scatterGlobalToLocal(solutionDotVec);

  // Set jacobian to zero.
  _jacobianLHS->zero();

  // Sum Jacobian contributions across integrators.
  const size_t numIntegrators = _integrators.size();
  for (size_t i=0; i < numIntegrators; ++i) {
    _integrators[i]->computeLHSJacobianImplicit(_jacobianLHS, _preconditionerLHS, t, dt, *_solution, *_solutionDot);
  } // for
  
  // Assemble jacobian.
  _jacobianLHS->assemble("final_assembly");

#if 0 // FIX THIS
  if (_customConstraintPCMat) {
    // Recalculate preconditioner.
    for (size_t i=0; i < numIntegrators; ++i) {
      _integrators[i]->computeLHSPreconditioner(&_customConstraintPCMat, _jacobianLHS, t, dt, *_solution, *_solutionDot);
    } // for

    MatAssemblyBegin(_customConstraintPCMat, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(_customConstraintPCMat, MAT_FINAL_ASSEMBLY);

#if 0 // debugging
    std::cout << "Preconditioner Matrix" << std::endl;
    MatView(_customConstraintPCMat, PETSC_VIEWER_STDOUT_WORLD);
#endif
  } // if
#endif

  PYLITH_METHOD_END;
} // computeLHSJacobianImplicit

// ----------------------------------------------------------------------
// Compute LHS Jacobian for F(t,s,\dot{s}) for explicit time stepping.
void
pylith::problems::Problem::computeLHSJacobianExplicit(const PylithReal t,
						      const PylithReal dt,
						      PetscVec solutionVec,
						      PetscVec solutionDotVec)
{ // computeLHSJacobianExplicit
  PYLITH_METHOD_BEGIN;

  assert(_jacobianLHS);
  assert(_solution);
  assert(_solutionDot);

  // :KLUDGE: Should add check to see if we need to compute Jacobian

  // Update PyLith view of the solution.
  _solution->scatterGlobalToLocal(solutionVec);
  _solutionDot->scatterGlobalToLocal(solutionDotVec);
  
  // Set jacobian to zero.
  _jacobianLHS->zero();

  // Sum Jacobian contributions across integrators.
  const size_t numIntegrators = _integrators.size();
  for (size_t i=0; i < numIntegrators; ++i) {
    _integrators[i]->computeLHSJacobianExplicit(_jacobianLHS, _preconditionerLHS, t, dt, *_solution, *_solutionDot);
  } // for
  
  // Assemble jacobian.
  _jacobianLHS->assemble("final_assembly");

  PYLITH_METHOD_END;
} // computeLHSJacobianExplicit

// End of file 
