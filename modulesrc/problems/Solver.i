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
 * @file modulesrc/problems/Solver.hh
 *
 * @brief Python interface to abstract base class Solver.
 */

namespace pylith {
  namespace problems {

    class Solver
    { // Solver

      // PUBLIC MEMBERS /////////////////////////////////////////////////
    public :

      /// Constructor.
      Solver(void);

      /// Destructor
      ~Solver(void);

      /** Initialize solver.
       *
       * @param fields Solution fields.
       * @param jacobian Jacobian of system.
       * @param formulation Formulation of system of equations.
       */
      virtual
      void
      initialize(const pylith::topology::SolutionFields& fields,
		 const pylith::topology::Jacobian& jacobian,
		 Formulation* const formulation);

    }; // Solver

  } // problems
} // pylith


// End of file 
