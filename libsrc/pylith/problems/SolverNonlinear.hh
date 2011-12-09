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
// Copyright (c) 2010-2011 University of California, Davis
//
// See COPYING for license information.
//
// ======================================================================
//

/**
 * @file libsrc/problems/SolverNonlinear.hh
 *
 * @brief Object for using PETSc scalable nonlinear equation solvers
 * (SNES).
 */

#if !defined(pylith_problems_solvernonlinear_hh)
#define pylith_problems_solvernonlinear_hh

// Include directives ---------------------------------------------------
#include "Solver.hh" // ISA Solver

#include "pylith/utils/petscfwd.h" // HASA PetscSNES
#include <petscmat.h> // USES MatStructure

// SolverNonlinear ---------------------------------------------------------
/** @brief Object for using PETSc scalable nonlinear equation solvers
 * (SNES).
 *
 * The PETSc nonlinear solvers provide an interface to Newton-based
 * methods for solving nonlinear equations.
 */
class pylith::problems::SolverNonlinear : public Solver
{ // SolverNonlinear
  friend class TestSolverNonlinear; // unit testing

// PUBLIC MEMBERS ///////////////////////////////////////////////////////
public :

  /// Constructor
  SolverNonlinear(void);

  /// Destructor
  ~SolverNonlinear(void);

  /// Deallocate PETSc and local data structures.
  void deallocate(void);
  
  /** Initialize solver.
   *
   * @param fields Solution fields.
   * @param jacobian Jacobian of system.
   * @param formulation Formulation of system of equations.
   */
  void
  initialize(const topology::SolutionFields& fields,
	     const topology::Jacobian& jacobian,
	     Formulation* const formulation);

  /** Solve the system.
   *
   * @param solveSoln Field to solve for.
   * @param jacobian Jacobian of the system.
   * @param residual Residual field.
   */
  void solve(topology::Field<topology::Mesh>* solveSoln,
	     topology::Jacobian* jacobian,
	     const topology::Field<topology::Mesh>& residual);

  /** Generic C interface for reformResidual for integration with
   * PETSc SNES solvers.
   *
   * @param snes PETSc scalable nonlinear equation solver.
   * @param tmpSolveSolnVec Temporary PETSc vector for solution.
   * @param tmpResidualVec Temporary PETSc vector for residual.
   * @param context ArgsResidual structure with arguments.
   * @returns PETSc error code.
   */
  static
  PetscErrorCode reformResidual(PetscSNES snes,
				PetscVec tmpSolveSolnVec,
				PetscVec tmpResidualVec,
				void* context);

  /** Generic C interface for reformJacobian for integration with
   * PETSc SNES solvers.
   *
   * @param snes PETSc scalable nonlinear equation solver.
   * @param tmpSolveSolnVec Temporary PETSc vector for solution.
   * @param jacobianMat PETSc sparse matrix for system Jacobian.
   * @param preconditionerMat PETSc sparse matrix for preconditioner.
   * @param Flag indicating layout of preconditioner matrix.
   * @param context ArgsJacobian structure with arguments.
   * @returns PETSc error code.
   */
  static
  PetscErrorCode reformJacobian(PetscSNES snes,
				PetscVec tmpSolveSolnVec,
				PetscMat* jacobianMat,
				PetscMat* preconditionerMat,
				MatStructure* preconditionerLayout,
				void* context);

  /** Generic C interface for customized PETSc line search.
   *
   * @param snes PETSc SNES solver.
   * @param lsctx Optional context for line search (not used here).
   * @param x Current iterate.
   * @param f Residual evaluated at x.
   * @param y Search direction.
   * @param w Work vector
   * @param f 2-norm of f.
   * @param xnorm Norm of x if known, otherwise 0.
   * @param g Residual evaluated at new iterate y.
   * @param w New iterate.
   * @param gnorm 2-norm of g.
   * @param ynorm 2-norm of search length.
   * @param PETSC_TRUE if line search succeeds; PETSC_FALSE on failure.
   * @returns PETSc error code.
   */
  static
  PetscErrorCode lineSearch(PetscSNES snes,
			    void *lsctx,
			    PetscVec x,
			    PetscVec f,
			    PetscVec y,
			    PetscReal fnorm,
			    PetscReal xnorm,
			    PetscVec g,
			    PetscVec w,
			    PetscReal *ynorm,
			    PetscReal *gnorm,
			    PetscBool *flag);

  /** Generic C interface for customized PETSc initial guess.
   *
   * @param snes PETSc SNES solver.
   * @param initialGuessVec PETSc vector for initial guess.
   * @param lsctx Optional context for line search (not used here).
   * @returns PETSc error code.
   */
  static
  PetscErrorCode initialGuess(PetscSNES snes,
			      PetscVec initialGuessVec,
			      void *lsctx);

// PRIVATE METHODS //////////////////////////////////////////////////////
private :

  /// Initialize logger.
  void _initializeLogger(void);

// PRIVATE MEMBERS //////////////////////////////////////////////////////
private :

  PetscSNES _snes; ///< PETSc SNES nonlinear solver.

// NOT IMPLEMENTED //////////////////////////////////////////////////////
private :

  SolverNonlinear(const SolverNonlinear&); ///< Not implemented
  const SolverNonlinear& operator=(const SolverNonlinear&); ///< Not implemented

}; // SolverNonlinear

#endif // pylith_problems_solvernonlinear_hh


// End of file 
