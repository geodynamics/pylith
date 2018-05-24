// -*- C++ -*-
//
// ----------------------------------------------------------------------
//
// Brad T. Aagaard, U.S. Geological Survey
// Charles A. Williams, GNS Science
// Matthew G. Knepley, University of Chicago
//
// This code was developed as part of the Computational Infrastructure
// for Geodynamics (http://geodynamics.org).
//
// Copyright (c) 2010-2017 University of California, Davis
//
// See COPYING for license information.
//
// ----------------------------------------------------------------------
//

/** @file modulesrc/faults/FaultCohesiveKin.i
 *
 * @brief Python interface to C++ FaultCohesiveKin object.
 */

namespace pylith {
  namespace faults {

      class FaultCohesiveKin : public FaultCohesive {

	  // PUBLIC METHODS /////////////////////////////////////////////////
      public :
	  
	  /// Default constructor.
	  FaultCohesiveKin(void);
	  
	  /// Destructor.
	  virtual
	  ~FaultCohesiveKin(void);
	  
	  /// Deallocate PETSc and local data structures.
	  void deallocate(void);
	  
	  /** Set kinematic earthquake sources.
	   *
	   * @param names Array of kinematic earthquake source names.
	   * @param numNames Number of earthquake sources.
	   * @param sources Array of kinematic earthquake sources.
	   * @param numSources Number of earthquake sources.
	   */
	  %apply(const char* const* string_list, const int list_len){
	      (const char* const* names,
	       const int numNames)
		  };
	  void eqsrcs(const char* const* names,
		      const int numNames,
		      pylith::faults::KinSrc** sources,
		      const int numSources);
	  %clear(const char* const* names, const int numNames);
	  
	  /** Initialize fault.
	   *
	   * Setup earthquake sources.
	   *
	   * @param[in] solution Solution field.
	   */
	  void initialize(const pylith::topology::Field& solution);
	  
	  /** Update auxiliary fields at beginning of time step.
	   *
	   * @param[in] t Current time.
	   * @param[in] dt Current time step.
	   */
	  void prestep(const double t,
		       const double dt);
	  
	  /** Compute RHS residual for G(t,s).
	   *
	   * @param[out] residual Field for residual.
	   * @param[in] t Current time.
	   * @param[in] dt Current time step.
	   * @param[in] solution Field with current trial solution.
	   */
	  void computeRHSResidual(pylith::topology::Field* residual,
				  const PylithReal t,
				  const PylithReal dt,
				  const pylith::topology::Field& solution);
	  
	  /** Compute RHS Jacobian and preconditioner for G(t,s).
	   *
	   * @param[out] jacobianMat PETSc Mat with Jacobian sparse matrix.
	   * @param[out] precondMat PETSc Mat with Jacobian preconditioning sparse matrix.
	   * @param[in] t Current time.
	   * @param[in] dt Current time step.
	   * @param[in] solution Field with current trial solution.
	   */
	  void computeRHSJacobian(PetscMat jacobianMat,
				  PetscMat preconMat,
				  const PylithReal t,
				  const PylithReal dt,
				  const pylith::topology::Field& solution);
	  
	  /** Compute LHS residual for F(t,s,\dot{s}).
	   *
	   * @param[out] residual Field for residual.
	   * @param[in] t Current time.
	   * @param[in] dt Current time step.
	   * @param[in] solution Field with current trial solution.
	   * @param[in] solutionDot Field with time derivative of current trial solution.
	   */
	  void computeLHSResidual(pylith::topology::Field* residual,
				  const PylithReal t,
				  const PylithReal dt,
				  const pylith::topology::Field& solution,
				  const pylith::topology::Field& solutionDot);
	  
	  /** Compute LHS Jacobian and preconditioner for F(t,s,\dot{s}) with implicit time-stepping.
	   *
	   * @param[out] jacobianMat PETSc Mat with Jacobian sparse matrix.
	   * @param[out] precondMat PETSc Mat with Jacobian preconditioning sparse matrix.
	   * @param[in] t Current time.
	   * @param[in] dt Current time step.
	   * @param[in] tshift Scale for time derivative.
	   * @param[in] solution Field with current trial solution.
	   * @param[in] solutionDot Field with time derivative of current trial solution.
	   */
	  void computeLHSJacobianImplicit(PetscMat jacobianMat,
					  PetscMat precondMat,
					  const PylithReal t,
					  const PylithReal dt,
					  const PylithReal tshift,
					  const pylith::topology::Field& solution,
					  const pylith::topology::Field& solutionDot);
	  
	  /** Compute inverse of lumped LHS Jacobian for F(t,s,\dot{s}) with explicit time-stepping.
	   *
	   * @param[out] jacobianInv Inverse of lumped Jacobian as a field.
	   * @param[in] t Current time.
	   * @param[in] dt Current time step.
	   * @param[in] tshift Scale for time derivative.
	   * @param[in] solution Field with current trial solution.
	   */
	  void computeLHSJacobianLumpedInv(pylith::topology::Field* jacobianInv,
					   const PylithReal t,
					   const PylithReal dt,
					   const PylithReal tshift,
					   const pylith::topology::Field& solution);
	  
	  
      }; // class FaultCohesiveKin
      
  } // faults
} // pylith


// End of file 
