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
 * @file modulesrc/problems/SolverLumped.hh
 *
 * @brief Python interface to abstract base class SolverLumped.
 */

namespace pylith {
  namespace problems {

    class SolverLumped
    { // SolverLumped

      // PUBLIC MEMBERS /////////////////////////////////////////////////
    public :

      /// Constructor.
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
      initialize(const pylith::topology::SolutionFields& fields,
		 const pylith::topology::Field& jacobian,
		 Formulation* const formulation);
      
      /** Solve the system.
       *
       * @param solution Solution field.
       * @param jacobian Jacobian of the system.
       * @param residual Residual field.
       */
      void solve(pylith::topology::Field* solution,
		 const pylith::topology::Field& jacobian,
		 const pylith::topology::Field& residual);
      
    }; // SolverLumped

  } // problems
} // pylith


// End of file 
