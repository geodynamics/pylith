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
 * @file pylith/problems/Solver.hh
 *
 * @brief Abstract C++ base class for using PETSc linear and nonlinear
 * solvers.
 */

#if !defined(pylith_problems_solver_hh)
#define pylith_problems_solver_hh

// Include directives ---------------------------------------------------
#include "problemsfwd.hh" // forward declarations

#include "pylith/topology/topologyfwd.hh" // USES SolutionFields
#include "pylith/utils/petscfwd.h" // USES PetscVec


// Solver ---------------------------------------------------------
class pylith::problems::Solver
{ // Integrator
  friend class TestSolver; // unit testing

// PUBLIC MEMBERS ///////////////////////////////////////////////////////
public :

  /// Constructor
  Solver(void);

  /// Destructor
  ~Solver(void);

  /** Initialize solver.
   *
   * @param fields Solution fields.
   */
  virtual
  void
  initialize(topology::SolutionFields* fields);

// PRIVATE MEMBERS //////////////////////////////////////////////////////
private :


// NOT IMPLEMENTED //////////////////////////////////////////////////////
private :

  Solver(const Solver&); ///< Not implemented
  const Solver& operator=(const Solver&); ///< Not implemented

}; // Solver

#endif // pylith_problems_solver_hh


// End of file 
