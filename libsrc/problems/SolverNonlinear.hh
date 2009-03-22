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

// SolverNonlinear ---------------------------------------------------------
class pylith::problems::SolverNonlinear
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
   * @param jacobian Jacobian of the system.
   * @param residualFn Function for reforming residual.
   * @param argsResidual Structure holding args for reforming residual.
   * @param jacobianFn Function for reforming residual.
   * @param argsJacobian Structure holding args for reforming residual.
   */
  void
  initialize(topology::SolutionFields* fields,
	     const topology::Jacobian& jacobian,
	     Formulation::ResidualFn* residualFn,
	     Formulation::ArgsResidual* argsResidual,
	     Formulation::JacobianFn* jacobianFn,
	     Formulation::ArgsJacobian* argsJacobian);

  /** Solve the system.
   *
   * @param solution Solution field.
   * @param jacobian Jacobian of the system.
   * @param residual Residual field.
   */
  void solve(topology::Field<topology::Mesh>* solution,
	     const topology::Jacobian& jacobian,
	     const topology::Field<topology::Mesh>& residual);

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
