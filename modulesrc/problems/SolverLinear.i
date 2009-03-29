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
 * @file modulesrc/problems/SolverLinear.hh
 *
 * @brief Python interface to abstract base class SolverLinear.
 */

namespace pylith {
  namespace problems {

    class SolverLinear
    { // SolverLinear

      // PUBLIC MEMBERS /////////////////////////////////////////////////
    public :

      /// Constructor.
      SolverLinear(void);

      /// Destructor
      ~SolverLinear(void);

      /// Deallocate PETSc and local data structures.
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
      initialize(const pylith::topology::SolutionFields& fields,
		 const pylith::topology::Jacobian& jacobian,
		 Formulation* const formulation);

      /** Solve the system.
       *
       * @param solution Solution field.
       * @param jacobian Jacobian of the system.
       * @param residual Residual field.
       */
      void solve(pylith::topology::Field<pylith::topology::Mesh>* solution,
		 const pylith::topology::Jacobian& jacobian,
		 const pylith::topology::Field<pylith::topology::Mesh>& residual);

    }; // SolverLinear

  } // problems
} // pylith


// End of file 
