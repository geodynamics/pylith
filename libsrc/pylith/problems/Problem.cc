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

#include "pylith/feassemble/Integrator.hh" // USES Integrator
#include "pylith/topology/Jacobian.hh" // USES Jacobian
#include "pylith/topology/SolutionFields.hh" // USES SolutionFields

#include "pylith/utils/error.hh" // USES PYLITH_CHECK_ERROR
#include <cassert> // USES assert()

// ----------------------------------------------------------------------
// Constructor
pylith::problems::Problem::Problem(void) :
  _solution(0),
  _jacobianRHS(0),
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

  _solution = 0; // :TODO: Use shared pointer.
  _jacobianRHS = 0; // :TODO: Use shared pointer.

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
pylith::problems::Problem::integrators(pylith::feassemble::Integrator* integratorArray[],
				       const int numIntegrators)
{ // integrators
  assert( (!integratorArray && 0 == numIntegrators) || (integratorArray && 0 < numIntegrators) );
  _integrators.resize(numIntegrators);
  for (int i=0; i < numIntegrators; ++i)
    _integrators[i] = integratorArray[i];
} // integrators
  
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
// Reform system residual.
void
pylith::problems::Problem::reformRHSResidual(const PetscVec* tmpResidualVec,
					     const PetscVec* tmpSolutionVec)
{ // reformRHSResidual
  PYLITH_METHOD_BEGIN;

  

  // Update section view of field.
  if (tmpSolutionVec) {
    topology::Field& solution = _fields->solution();
    solution.scatterGlobalToLocal(*tmpSolutionVec);
  } // if

  // Update rate fields (must be consistent with current solution).
  calcRateFields();  

  // Set residual to zero.
  topology::Field& residual = _fields->get("residual");
  residual.zeroAll();

  // Add in contributions that require assembly.
  const int numIntegrators = _integrators.size();
  assert(numIntegrators > 0); // must have at least 1 integrator
  for (int i=0; i < numIntegrators; ++i) {
    _integrators[i]->timeStep(_dt);
    _integrators[i]->integrateResidual(residual, _t, _fields);
  } // for

  // Assemble residual.
  residual.complete();

  // Update PETSc view of residual
  if (tmpResidualVec)
    residual.scatterLocalToGlobal(*tmpResidualVec);
  else
    residual.scatterLocalToGlobal();

  // TODO: Move this to SolverLinear 
  if (tmpResidualVec)
    VecScale(*tmpResidualVec, -1.0);

  PYLITH_METHOD_END;
} // reformRHSResidual

// ----------------------------------------------------------------------
// Reform system Jacobian.
void
pylith::problems::Problem::reformRHSJacobian(const PetscVec* tmpSolutionVec)
{ // reformRHSJacobian
  PYLITH_METHOD_BEGIN;

  assert(_jacobianRHS);
  assert(_solution);

  // Update section view of field.
  if (tmpSolutionVec) {
    _solution->scatterGlobalToLocal(*tmpSolutionVec);
  } // if

  // Set jacobian to zero.
  _jacobianRHS->zero();

  // Add in contributions that require assembly.
  const int numIntegrators = _integrators.size();
  for (int i=0; i < numIntegrators; ++i) {
    _integrators[i]->integrateRHSJacobian(_jacobianRHS, _t, _fields);
  } // for
  
  // Assemble jacobian.
  _jacobianRHS->assemble("final_assembly");

  if (_customConstraintPCMat) {
    // Recalculate preconditioner.
    for (int i=0; i < numIntegrators; ++i) {
      _integrators[i]->integrateRHSPreconditioner(&_customConstraintPCMat, _jacobianRHS, *_solution);
    } // for

    MatAssemblyBegin(_customConstraintPCMat, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(_customConstraintPCMat, MAT_FINAL_ASSEMBLY);

#if 0 // debugging
    std::cout << "Preconditioner Matrix" << std::endl;
    MatView(_customConstraintPCMat, PETSC_VIEWER_STDOUT_WORLD);
#endif
  } // if

  PYLITH_METHOD_END;
} // reformRHSJacobian

// End of file 
