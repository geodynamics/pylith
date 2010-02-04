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
 * @file libsrc/problems/SolverLumped.hh
 *
 * @brief Object for using PETSc scalable linear equation solvers (KSP).
 */

#if !defined(pylith_problems_solverlumped_hh)
#define pylith_problems_solverlumped_hh

// Include directives ---------------------------------------------------
#include "Solver.hh" // ISA Solver

#include "pylith/utils/petscfwd.h" // HASA PetscKSP

// SolverLumped ---------------------------------------------------------
/** @brief Object for using simple solver to solver system with lumped Jacobian.
 */

class pylith::problems::SolverLumped : Solver
{ // SolverLumped
  friend class TestSolverLumped; // unit testing

// PUBLIC MEMBERS ///////////////////////////////////////////////////////
public :

  /// Constructor
  SolverLumped(void);

  /// Destructor
  ~SolverLumped(void);

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
	     const topology::Field<topology::Mesh>& jacobian,
	     Formulation* const formulation);

  /** Solve the system.
   *
   * @param solution Solution field.
   * @param jacobian Jacobian of the system.
   * @param residual Residual field.
   */
  void solve(topology::Field<topology::Mesh>* solution,
	     const topology::Field<topology::Mesh>& jacobian,
	     const topology::Field<topology::Mesh>& residual);

// PRIVATE METHODS //////////////////////////////////////////////////////
private :

  /// Initialize logger.
  void _initializeLogger(void);

// NOT IMPLEMENTED //////////////////////////////////////////////////////
private :

  SolverLumped(const SolverLumped&); ///< Not implemented
  const SolverLumped& operator=(const SolverLumped&); ///< Not implemented

}; // SolverLumped

#endif // pylith_problems_solverlumped_hh


// End of file 
