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
 * @file modulesrc/problems/Problem.hh
 *
 * @brief Python interface to C++ Problem.
 */

namespace pylith {
  namespace problems {

    class Problem
    { // Problem

      // PUBLIC ENUM ////////////////////////////////////////////////////
    public :
      
      enum SolverTypeEnum {
	LINEAR, // Linear solver.
	NONLINEAR, // Nonlinear solver.
      }; // SolverType
      
      
      // PUBLIC MEMBERS /////////////////////////////////////////////////
    public :
      
      /// Constructor
      Problem(void);
      
      /// Destructor
      virtual
      ~Problem(void);
      
      /// Deallocate PETSc and local data structures.
      void deallocate(void);
      
      /** Set solver type.
       *
       * @param[in] value Solver type.
       */
      void solverType(const SolverTypeEnum value);
      
      /** Get solver type.
       *
       * @returns Solver type.
       */
      SolverTypeEnum solverType(void) const;

      /** Set handles to integrators.
       *
       * @param[in] integratorArray Array of integrators.
       * @param[in] numIntegrators Number of integrators.
       */
      void integrators(pylith::feassemble::IntegratorPointwise* integratorArray[] ,
		       const int numIntegrators);
      
      /** Set handles to constraints.
       *
       * @param[in] constraintArray Array of constraints.
       * @param[in] numContraints Number of constraints.
       */
      void constraints(pylith::feassemble::Constraint* constraintArray[] ,
		       const int numConstraints);
      
      /** Initialize.
       *
       */
      virtual
      void initialize(void);
      
      /** Compute RHS residual, G(t,s).
       *
       * @param[out] residualVec PETSc Vec for residual.
       * @param[in] t Current time.
       * @param[in] dt Current time step.
       * @param[in] solutionVec PETSc Vec with current trial solution.
       */
      void computeRHSResidual(PetscVec residualVec,
			      const PetscReal t,
			      const PetscReal dt,
			      PetscVec solutionVec);
      
      /* Compute RHS Jacobian for G(t,s).
       *
       * @param[out] jacobianMat PETSc Mat for Jacobian.
       * @param[out] precondMat PETSc Mat for preconditioner for Jacobian.
       * @param[in] t Current time.
       * @param[in] dt Current time step.
       * @param[in] solutionVec PETSc Vec with current trial solution.
       */
      void computeRHSJacobian(PetscMat jacobianMat,
			      PetscMat precondMat,
			      const PylithReal t,
			      const PylithReal dt,
			      PetscVec solutionVec);
      
      /** Compute LHS residual, F(t,s,\dot{s}).
       *
       * @param[out] residualVec PETSc Vec for residual.
       * @param[in] t Current time.
       * @param[in] dt Current time step.
       * @param[in] solutionVec PETSc Vec with current trial solution.
       * @param[in] solutionDotVec PETSc Vec with time derivative of current trial solution.
       */
      void computeLHSResidual(PetscVec residualVec,
			      const PetscReal t,
			      const PetscReal dt,
			      PetscVec solutionVec,
			      PetscVec solutionDotVec);
      
      /* Compute LHS Jacobian for F(t,s,\dot{s}) for implicit time stepping.
       *
       * @param[out] jacobianMat PETSc Mat for Jacobian.
       * @param[out] precondMat PETSc Mat for preconditioner for Jacobian.
       * @param[in] t Current time.
       * @param[in] dt Current time step.
       * @param[in] tshift Scale for time derivative.
       * @param[in] solutionVec PETSc Vec with current trial solution.
       * @param[in] solutionDotVec PETSc Vec with time derivative of current trial solution.
       */
      void computeLHSJacobianImplicit(PetscMat jacobianMat,
				      PetscMat precondMat,
				      const PylithReal t,
				      const PylithReal dt,
				      const PylithReal tshift,
				      PetscVec solutionVec,
				      PetscVec solutionDotVec);
      
      /* Compute inverse of lumped LHS Jacobian for F(t,s,\dot{s}) for explicit time stepping.
       *
       * @param[in] t Current time.
       * @param[in] dt Current time step.
       * @param[in] solutionVec PETSc Vec with current trial solution.
       */
      void computeLHSJacobianLumpedInv(const PylithReal t,
				       const PylithReal dt,
				       PetscVec solutionVec);
      
      
    }; // Problem
    
  } // problems
} // pylith


// End of file
