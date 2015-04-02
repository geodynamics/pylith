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
      initialize(const pylith::topology::SolutionFields& fields,
		 const pylith::topology::Jacobian& jacobian,
		 Formulation* const formulation);

    }; // Solver

  } // problems
} // pylith


// End of file 
