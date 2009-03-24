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
 * @file modulesrc/problems/SolverNonlinear.hh
 *
 * @brief Python interface to abstract base class SolverNonlinear.
 */

namespace pylith {
  namespace problems {

    class SolverNonlinear
    { // SolverNonlinear

      // PUBLIC MEMBERS /////////////////////////////////////////////////
    public :

      /// Constructor.
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

    }; // SolverNonlinear

  } // problems
} // pylith


// End of file 
