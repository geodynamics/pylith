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
 * @file libsrc/problems/SolverLinear.hh
 *
 * @brief Object for using PETSc scalable linear equation solvers (KSP).
 *
 * The PETSc linear KSP solvers provide an interface to Krylov subspace
 * (KS) iterative methods and preconditioners (P).
 */

#if !defined(pylith_problems_solverlinear_hh)
#define pylith_problems_solverlinear_hh

// Include directives ---------------------------------------------------
#include "Solver.hh" // ISA Solver

#include "pylith/utils/petscfwd.h" // HASA PetscKSP

// SolverLinear ---------------------------------------------------------
class pylith::problems::SolverLinear : Solver
{ // SolverLinear
  friend class TestSolverLinear; // unit testing

// PUBLIC MEMBERS ///////////////////////////////////////////////////////
public :

  /// Constructor
  SolverLinear(void);

  /// Destructor
  ~SolverLinear(void);

  /// Deallocate data structures.
  void deallocate(void);
  
  /** Set initial guess zero flag.
   *
   * @param value True means use zero as initial guess, false means
   * use previous solution as initial guess.
   */
  void initialGuessZero(const bool value);

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

// PRIVATE MEMBERS //////////////////////////////////////////////////////
private :

  PetscKSP _ksp; ///< PETSc KSP linear solver.

// NOT IMPLEMENTED //////////////////////////////////////////////////////
private :

  SolverLinear(const SolverLinear&); ///< Not implemented
  const SolverLinear& operator=(const SolverLinear&); ///< Not implemented

}; // SolverLinear

#endif // pylith_problems_solverlinear_hh


// End of file 
