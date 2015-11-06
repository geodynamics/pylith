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
 * @file libsrc/problems/TimeDependent.hh
 *
 * @brief Object for time dependent problem.
 */

#if !defined(pylith_problems_timedependent_hh)
#define pylith_problems_timedependent_hh

// Include directives ---------------------------------------------------
#include "Problem.hh" // ISA Problem

// TimeDependent ---------------------------------------------------------
/** @brief Object for time dependent problem.
 *
 * TimeDependent uses the PETSc TS object for time stepping.
 */

class pylith::problems::TimeDependent : public Problem
{ // TimeDependent
  friend class TestTimeDependent; // unit testing

// PUBLIC ENUM //////////////////////////////////////////////////////////
public :

  enum ProblemTypeEnum {
    LINEAR, // Linear solver.
    NONLINEAR, // Nonlinear solver.
  }; // ProblemType

  enum FormulationTypeEnum {
    IMPLICIT, // Implicit time stepping.
    EXPLICIT, // Explicit time stepping.
  }; // FormulationTypeEnum

// PUBLIC MEMBERS ///////////////////////////////////////////////////////
public :

  /// Constructor
  TimeDependent(void);

  /// Destructor
  ~TimeDependent(void);

  /// Deallocate PETSc and local data structures.
  void deallocate(void);

  /** Set problem type.
   *
   * @param value Problem type.
   */
  void problemType(const ProblemTypeEnum value);

  /** Get problem type.
   *
   * @returns Problem type.
   */
  ProblemTypeEnum problemType(void) const;

  /** Set start time for problem.
   *
   * @param value Start time (nondimensional).
   */
  void startTime(const PetscReal value);

  /** Get start time for problem.
   *
   * @returns Start time (nondimensional).
   */
  PetscReal startTime(void) const;

  /** Set total time for problem.
   *
   * @param value Total time (nondimensional).
   */
  void totalTime(const PetscReal value);

  /** Get total time for problem.
   *
   * @returns Total time (nondimensional).
   */
  PetscReal totalTime(void) const;

  /** Set maximum number of time steps.
   *
   * @param value Maximum number of time steps.
   */
  void maxTimeSteps(const PetscInt value);

  /** Get maximum number of time steps.
   *
   * @returns Maximum number of time steps.
   */
  PetscInt maxTimeSteps(void) const;

  /** Set initial time step for problem.
   *
   * @param value Initial time step (nondimensional).
   */
  void dtInitial(const PetscReal value);

  /** Get initial time step for problem.
   *
   * @returns Initial time step (nondimensional).
   */
  PetscReal dtInitial(void) const;

  /** Initialize.
   */
  void initialize(pylith::topology::Field* solution,
		  pylith::topology::Jacobian* jacobianRHS);
  
  /** Solve time dependent problem.
   */
  void solve(void);
  
  /** Perform operations before advancing solution one time step.
   *
   * Set constraints, etc.
   */
  void prestep(void);

  /** Perform Perform operations after advancing solution one time step.
   *
   * Update state variables, output.
   */
  void poststep(void);

  /** Callback static method for computing residual for RHS, G(t,s).
   *
   * @param ts PETSc time stepper.
   * @param t Current time.
   * @param solutionVec PetscVec for solution.
   * @param residualvec PetscVec for residual.
   * @param context User context (TimeDependent).
   */
  static
  PetscErrorCode computeRHSResidual(PetscTS ts,
				    PetscReal t,
				    PetscVec solutionVec,
				    PetscVec residualVec,
				    void* context);
  
  /* Callback static method for computing Jacobian for RHS, Jacobian of G(t,s).
   *
   * @param ts PETSc time stepper.
   * @param t Current time.
   * @param solution PetscVec for solution.
   * @param jacobianMat Jacobian matrix.
   * @param precondMat Preconditioner matrix.
   * @param context User context (TimeDependent).
   */
  static
  PetscErrorCode computeRHSJacobian(PetscTS ts,
				    PetscReal t,
				    PetscVec solutionVec,
				    PetscMat jacobianMat,
				    PetscMat precondMat,
				    void* context);

  /** Callback static method for computing residual for LHS, F(t,s,\dot{s}).
   *
   * @param ts PETSc time stepper.
   * @param t Current time.
   * @param solutionVec PetscVec for solution.
   * @param solutionDotVec PetscVec for time derivative of solution.
   * @param residualvec PetscVec for residual.
   * @param context User context (TimeDependent).
   */
  static
  PetscErrorCode computeLHSResidual(PetscTS ts,
				    PetscReal t,
				    PetscVec solutionVec,
				    PetscVec solutionDotVec,
				    PetscVec residualVec,
				    void* context);
  
  /* Callback static method for computing Jacobian for LHS, Jacobian of F(t,s,\dot{s}).
   *
   * @param ts PETSc time stepper.
   * @param t Current time.
   * @param solution PetscVec for solution.
   * @param jacobianMat Jacobian matrix.
   * @param precondMat Preconditioner matrix.
   * @param context User context (TimeDependent).
   */
  static
  PetscErrorCode computeLHSJacobian(PetscTS ts,
				    PetscReal t,
				    PetscVec solutionVec,
				    PetscVec solutionDotVec,
				    PetscReal tshift,
				    PetscMat jacobianMat,
				    PetscMat precondMat,
				    void* context);

  /** Callback static method for operations before advancing solution one time step.
   */
  static
  PetscErrorCode prestep(PetscTS ts);

  /** Callback static method for operations after advancing solution one time step.
   */
  static
  PetscErrorCode poststep(PetscTS ts);
  
// PRIVATE MEMBERS //////////////////////////////////////////////////////
private :

  PetscReal _startTime; ///< Starting time.
  PetscReal _dtInitial; ///< Initial time step.
  PetscReal _totalTime; ///< Total time (duration) of problem.
  PetscInt _maxTimeSteps; ///< Maximum number of time steps for problem.
  PetscTS _ts; ///< PETSc time stepper.
  ProblemTypeEnum _problemType; ///< Problem (solver) type.
  FormulationTypeEnum _formulationType; ///< Type of time stepping.

// NOT IMPLEMENTED //////////////////////////////////////////////////////
private :

  TimeDependent(const TimeDependent&); ///< Not implemented
  const TimeDependent& operator=(const TimeDependent&); ///< Not implemented

}; // TimeDependent

#endif // pylith_problems_timedependent_hh


// End of file 
