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
 * @file libsrc/problems/Problem.hh
 *
 * @brief C++ object that manages formulating the equations.
 */

#if !defined(pylith_problems_problem_hh)
#define pylith_problems_problem_hh

// Include directives ---------------------------------------------------
#include "problemsfwd.hh" // forward declarations

#include "pylith/feassemble/feassemblefwd.hh" // USES Integrator
#include "pylith/topology/topologyfwd.hh" // USES Mesh, Field

#include "pylith/utils/petscfwd.h" // USES PetscVec, PetscMat

#include "pylith/utils/array.hh" // HASA std::vector


// Problem ----------------------------------------------------------
/** Reform the Jacobian and residual for the problem.
 *
 * We cast the problem in terms of F(t,s,\dot{s}) = G(t,s), s(t0) = s0.
 *
 * In PETSc time stepping (TS) notation, G is the RHS, and F is the I
 * function (which we call the LHS).
 *
 */
class pylith::problems::Problem
{ // Problem
  friend class TestProblem; // unit testing

// PUBLIC MEMBERS ///////////////////////////////////////////////////////
public :

  /// Constructor
  Problem(void);

  /// Destructor
  ~Problem(void);

  /// Deallocate PETSc and local data structures.
  void deallocate(void);

  /** Set handles to integrators.
   *
   * @param integratorArray Array of integrators.
   * @param numIntegrators Number of integrators.
   */
  void integrators(feassemble::IntegratorPointwise* integratorArray[] ,
		   const int numIntegrators);
  
  /** Set handles to constraints.
   *
   * @param constraintArray Array of constraints.
   * @param numContraints Number of constraints.
   */
  void constraints(feassemble::Constraint* constraintArray[] ,
		   const int numConstraints);
  
  /** Set handle to preconditioner.
   *
   * @param pc PETSc preconditioner.
   */
  void customPCMatrix(PetscMat& mat);

  /** Initialize.
   *
   */
  virtual
  void initialize(void);

  /** Compute RHS residual, G(t,s).
   *
   * @param[in] t Current time.
   * @param[in] dt Current time step.
   * @param[in] solutionVec PETSc Vec with current trial solution.
   * @param[out] residualVec PETSc Vec for residual.
   */
  void computeRHSResidual(const PetscReal t,
			  const PetscReal dt,
			  PetscVec solutionVec,
			  PetscVec residualVec);
  
  /* Compute RHS Jacobian for G(t,s).
   *
   * @param[in] t Current time.
   * @param[in] dt Current time step.
   * @param[in] solutionVec PETSc Vec with current trial solution.
   * @param[out] jacobianMat PETSc Mat for Jacobian. 
   * @param[out] precondMat PETSc Mat for preconditioner for Jacobian. 
   */
  void computeRHSJacobian(const PylithReal t,
			  const PylithReal dt,
			  PetscVec solutionVec,
			  PetscMat jacobianMat,
			  PetscMat precondMat);

  /** Compute LHS residual, F(t,s,\dot{s}).
   *
   * @param t Current time.
   * @param dt Current time step.
   * @param solutionVec PETSc Vec with current trial solution.
   * @param[in] solutionDotVec PETSc Vec with time derivative of current trial solution.
   * @param residualVec PETSc Vec for residual.
   */
  void computeLHSResidual(const PetscReal t,
			  const PetscReal dt,
			  PetscVec solutionVec,
			  PetscVec solutionDotVec,
			  PetscVec residualVec);
  
  /* Compute LHS Jacobian for F(t,s,\dot{s}) for implicit time stepping.
   *
   * @param[in] t Current time.
   * @param[in] dt Current time step.
   * @param[in] tshift Scale for time derivative.
   * @param[in] solutionVec PETSc Vec with current trial solution.
   * @param[in] solutionDotVec PETSc Vec with time derivative of current trial solution.
   * @param[out] jacobianMat PETSc Mat for Jacobian. 
   * @param[out] precondMat PETSc Mat for preconditioner for Jacobian. 
   */
  void computeLHSJacobianImplicit(const PylithReal t,
				  const PylithReal dt,
				  const PylithReal tshift,
				  PetscVec solutionVec,
				  PetscVec solutionDotVec,
				  PetscMat jacobianMat,
				  PetscMat precondMat);

  /* Compute LHS Jacobian for F(t,s,\dot{s}) for explicit time stepping.
   *
   * @param[in] t Current time.
   * @param[in] dt Current time step.
   * @param[in] solutionVec PETSc Vec with current trial solution.
   */
  void computeLHSJacobianExplicit(const PylithReal t,
				  const PylithReal dt,
				  PetscVec solutionVec);

// PROTECTED MEMBERS ////////////////////////////////////////////////////
protected :

  pylith::topology::Field* _solution; ///< Handle to solution field.
  pylith::topology::Field* _residualRHS; ///< Handle to residual field for RHS, G(t,s).
  pylith::topology::Field* _residualLHS; ///< Handle to residual field for LHS, F(t,s,\dot{s}).
  pylith::topology::Jacobian* _jacobianRHS; ///< Handle to Jacobian for RHS, G(t,s).
  pylith::topology::Jacobian* _jacobianLHS; ///< Handle to Jacobian for LHS, F(t,s,\dot{s}).
  pylith::topology::Jacobian* _preconditionerRHS; ///< Handle to Jacobian preconditioner for RHS, G(t,s).
  pylith::topology::Jacobian* _preconditionerLHS; ///< Handle to Jacobian preconditioner for LHS, F(t,s,\dot{s}).

  std::vector<pylith::feassemble::IntegratorPointwise*> _integrators; ///< Array of integrators.
  std::vector<pylith::feassemble::Constraint*> _constraints; ///< Array of constraints.

  PetscMat _customConstraintPCMat; ///< Custom PETSc preconditioning matrix for constraints.
  bool _useCustomConstraintPC; ///< True if using custom preconditioner for Lagrange constraints.

// NOT IMPLEMENTED //////////////////////////////////////////////////////
private :

  Problem(const Problem&); ///< Not implemented
  const Problem& operator=(const Problem&); ///< Not implemented

}; // Problem

#endif // pylith_problems_problem_hh


// End of file 
