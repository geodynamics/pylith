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
// Copyright (c) 2010-2015 University of California, Davis
//
// See COPYING for license information.
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

      /// Deallocate PETSc and local data structures.
      void deallocate(void);

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
      void solve(pylith::topology::Field* solution,
		 pylith::topology::Jacobian* jacobian,
		 const pylith::topology::Field& residual);

    }; // SolverNonlinear

  } // problems
} // pylith


// End of file 
