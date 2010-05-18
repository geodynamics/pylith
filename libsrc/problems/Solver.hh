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
 * @file libsrc/problems/Solver.hh
 *
 * @brief Abstract C++ base class for using PETSc linear and nonlinear
 * solvers.
 */

#if !defined(pylith_problems_solver_hh)
#define pylith_problems_solver_hh

// Include directives ---------------------------------------------------
#include "problemsfwd.hh" // forward declarations

#include "pylith/topology/topologyfwd.hh" // USES SolutionFields
#include "pylith/utils/utilsfwd.hh" // USES EventLogger
#include "pylith/utils/petscfwd.h" // USES PetscPC

// Solver ---------------------------------------------------------
/** @brief Abstract C++ base class for using PETSc linear and
 * nonlinear solvers.
 */
class pylith::problems::Solver
{ // Solver
  friend class TestSolver; // unit testing

// PUBLIC MEMBERS ///////////////////////////////////////////////////////
public :

  /// Constructor.
  Solver(void);

  /// Destructor
  virtual
  ~Solver(void);

  /// Deallocate PETSc and local data structures.
  virtual
  void deallocate(void);
  
  /** Initialize solver.
   *
   * @param fields Solution fields.
   * @param jacobian Jacobian of system.
   * @param formulation Formulation of system of equations.
   */
  virtual
  void
  initialize(const topology::SolutionFields& fields,
	     const topology::Jacobian& jacobian,
	     Formulation* const formulation);

// PROTECTED METHODS ////////////////////////////////////////////////////
protected :

  /** Setup preconditioner for preconditioning with split fields.8
   *
   * @param pc PETSc preconditioner.
   * @param formulation Formulation of system of equations.
   * @param fields Solution fields.
   */
  void
  _setupFieldSplit(PetscPC* const pc,
		   Formulation* const formulation,
		   const topology::SolutionFields& fields);

// PROTECTED MEMBERS ////////////////////////////////////////////////////
protected :

  Formulation* _formulation; ///< Handle to formulation for system of eqns.
  utils::EventLogger* _logger; ///< Event logger.

// NOT IMPLEMENTED //////////////////////////////////////////////////////
private :

  Solver(const Solver&); ///< Not implemented
  const Solver& operator=(const Solver&); ///< Not implemented

}; // Solver

#endif // pylith_problems_solver_hh


// End of file 
