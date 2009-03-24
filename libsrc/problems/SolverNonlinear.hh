// -*- C++ -*-
//
// ======================================================================
//
//                           Brad T. Aagaard
//                        U.S. Geological Survey
//
// {LicenseText}
//
// ======================================================================
//

/**
 * @file pylith/problems/SolverNonlinear.hh
 *
 * @brief Object for using PETSc scalable nonlinear equation solvers
 * (SNES).
 *
 * The PETSc nonlinear solvers provide an interface to Newton-based
 * methods for solving nonlinear equations.
 */

#if !defined(pylith_problems_solvernonlinear_hh)
#define pylith_problems_solvernonlinear_hh

// Include directives ---------------------------------------------------
#include "Solver.hh" // ISA Solver

#include "pylith/utils/petscfwd.h" // HASA PetscSNES
#include <petscmat.h> // USES MatStructure

// SolverNonlinear ---------------------------------------------------------
class pylith::problems::SolverNonlinear : public Solver
{ // Integrator
  friend class TestSolverNonlinear; // unit testing

// PUBLIC MEMBERS ///////////////////////////////////////////////////////
public :

  /// Constructor
  SolverNonlinear(void);

  /// Destructor
  ~SolverNonlinear(void);

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
   * @param solution Solution field.
   * @param jacobian Jacobian of the system.
   * @param residual Residual field.
   */
  void solve(topology::Field<topology::Mesh>* solution,
	     const topology::Jacobian& jacobian,
	     const topology::Field<topology::Mesh>& residual);

  /** Generic C interface for reformResidual for integration with
   * PETSc SNES solvers.
   *
   * @param snes PETSc scalable nonlinear equation solver.
   * @param solutionVec PETSc vector for solution.
   * @param residualVec PETSc vector for residual.
   * @param context ArgsResidual structure with arguments.
   * @returns PETSc error code.
   */
  static
  PetscErrorCode reformResidual(PetscSNES snes,
				PetscVec solutionVec,
				PetscVec residualVec,
				void* context);

  /** Generic C interface for reformJacobian for integration with
   * PETSc SNES solvers.
   *
   * @param snes PETSc scalable nonlinear equation solver.
   * @param solutionVec PETSc vector for solution.
   * @param jacobianMat PETSc sparse matrix for system Jacobian.
   * @param preconditionerMat PETSc sparse matrix for preconditioner.
   * @param Flag indicating layout of preconditioner matrix.
   * @param context ArgsJacobian structure with arguments.
   * @returns PETSc error code.
   */
  static
  PetscErrorCode reformJacobian(PetscSNES snes,
				PetscVec solutionVec,
				PetscMat* jacobianMat,
				PetscMat* preconditionerMat,
				MatStructure* preconditionerLayout,
				void* context);

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
