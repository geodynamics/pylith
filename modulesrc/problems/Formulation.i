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
// Copyright (c) 2010-2016 University of California, Davis
//
// See COPYING for license information.
//
// ======================================================================
//

/**
 * @file modulesrc/problems/Formulation.hh
 *
 * @brief Python interface to C++ Formulation.
 */

namespace pylith {
  namespace problems {

    class Formulation
    { // Formulation
      
      // PUBLIC MEMBERS /////////////////////////////////////////////////
    public :

      /// Constructor
      Formulation(void);

      /// Destructor
      ~Formulation(void);

      /// Deallocate PETSc and local data structures.
      void deallocate(void);
  
      /** Set flag for splitting fields.
       *
       * @param flag True if splitting fields, false otherwise.
       */
      void splitFields(const bool flag);

      /** Get flag for splitting fields.
       *
       * @returns flag True if splitting fields, false otherwise.
       */
      bool splitFields(void) const;

      /** Set flag for using custom Lagrange constraint preconditioner.
       *
       * @param flag True if using custom constraint precondition,
       * false otherwise.
       */
      void useCustomConstraintPC(const bool flag);

      /** Get flag indicating use of custom Lagrange constraint conditioner.
       *
       * @returns flag True if using custom constraint preconditioner, 
       * false otherwise.
       */
      bool useCustomConstraintPC(void) const;

      /** Get solution fields.
       *
       * @returns solution fields.
       */
      const pylith::topology::SolutionFields& fields(void) const;

      /** Get flag indicating whether Jacobian is symmetric.
       *
       * @returns True if Jacobian is symmetric, otherwise false.
       */
      bool isJacobianSymmetric(void) const;
  
      /** Set handles to integrators.
       *
       * @param integratorArray Array of integrators.
       * @param numIntegrators Number of integrators.
       */
      void integrators(pylith::feassemble::Integrator* integratorArray[],
		       const int numIntegrators);
      
      /** Update handles and parameters for reforming the Jacobian and
       *  residual.
       *
       * @param jacobian Handle to sparse matrix for Jacobian of system.
       * @param fields Handle to solution fields.
       * @param t Current time (nondimensional).
       * @param dt Time step (nondimension).
       */
      void updateSettings(pylith::topology::Jacobian* jacobian,
			  pylith::topology::SolutionFields* fields,
			  const PylithScalar t,
			  const PylithScalar dt);
      
      /** Update handles and parameters for reforming the Jacobian and
       *  residual.
       *
       * @param jacobian Handle to diagonal matrix (as Field) for
       * system Jacobian.
       * @param fields Handle to solution fields.
       * @param t Current time (nondimensional).
       * @param dt Time step (nondimension).
       */
      void updateSettings(pylith::topology::Field* jacobian,
			  pylith::topology::SolutionFields* fields,
			  const PylithScalar t,
			  const PylithScalar dt);

      /** Reform system residual.
       *
       * @param tmpResidualVec Temporary PETSc vector for residual.
       * @param tmpSolveSolnVec Temporary PETSc vector for solution.
       */
      void reformResidual(const PetscVec* tmpResidualVec =0,
			  const PetscVec* tmpSolveSolnVec =0);
      
      /* Reform system Jacobian.
       *
       * @param tmpSolveSolnVec Temporary PETSc vector for solution.
       */
      void reformJacobian(const PetscVec* tmpSolveSolnVec =0);
      
      /* Reform system Jacobian.
       */
      void reformJacobianLumped(void);

      /** Constrain solution space.
       *
       * @param tmpSolutionVec Temporary PETSc vector for solution.
       */
      void constrainSolnSpace(const PetscVec* tmpSolutionVec);

      /// Compute rate fields (velocity and/or acceleration) at time t.
      virtual
      void calcRateFields(void) = 0;

    }; // Formulation

  } // problems
} // pylith


// End of file 
