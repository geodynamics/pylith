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
 * @file libsrc/feassemble/ElasticityImplicitLgDeform.hh
 *
 * @brief Implicit time integration of quasi-static elasticity equation
 * with large rigid body motion and small strains using finite-elements.
 *
 * Note: This object operates on a single finite-element family, which
 * is defined by the quadrature and a database of material property
 * parameters.
 *
 * Computes contributions to terms A and r in
 *
 * A(t+dt) du(t) = b(t+dt, u(t), u(t-dt)) - A(t+dt) u(t)
 *
 * r = b(t+dt) - A(t+dt) (u(t) + du(t))
 *
 * where A(t) is a sparse matrix or vector, u(t+dt) is the field we
 * want to compute at time t+dt, b(t+dt) is a vector that depends on the
 * field at time t+dt and t, and du is zero at unknown DOF and set to
 * the constrained values at known DOF.
 *
 * Contributions to the RHS (b) include body forces, which are either
 * independent of u (small strain case) or are computed based on the
 * displacements at time t. The RHS also includes the internal force
 * vector, which is either constant for a time step (small strain,
 * elastic rheology) or changes with each iteration (large strain or
 * non-elastic rheology). The internal force vector is subtracted from the
 * existing force vector to get the residual. This object also computes
 * the entire stiffness matrix (A).
 *
 * See governing equations section of user manual for more
 * information.
 */

#if !defined(pylith_feassemble_elasticityimplicitlgdeform_hh)
#define pylith_feassemble_elasticityimplicitlgdeform_hh

// Include directives ---------------------------------------------------
#include "IntegratorElasticityLgDeform.hh" // ISA IntegratorElasticityLgDeform

// ElasticityImplicitLgDeform -------------------------------------------
class pylith::feassemble::ElasticityImplicitLgDeform : 
  public IntegratorElasticityLgDeform
{ // ElasticityImplicitLgDeform
  friend class TestElasticityImplicitLgDeform; // unit testing

// PUBLIC MEMBERS ///////////////////////////////////////////////////////
public :

  /// Constructor
  ElasticityImplicitLgDeform(void);

  /// Destructor
  ~ElasticityImplicitLgDeform(void);

  /// Deallocate PETSc and local data structures.
  void deallocate(void);
  
  /** Set time step for advancing from time t to time t+dt.
   *
   * @param dt Time step
   */
  void timeStep(const PylithScalar dt);

  /** Get stable time step for advancing from time t to time t+dt.
   *
   * Default is current time step.
   *
   * @param mesh Finite-element mesh.
   * @returns Time step
   */
  PylithScalar stableTimeStep(const topology::Mesh& mesh) const;

  /** Integrate residual part of RHS for 3-D finite elements.
   * Includes gravity and element internal force contribution.
   *
   * We assume that the effects of boundary conditions are already
   * included in the residual (tractions, concentrated nodal forces,
   * and contributions to internal force vector due to
   * displacement/velocity BC).  This routine computes the additional
   * external loads due to body forces plus the
   * element internal forces for the current stress state.
   *
   * @param residual Field containing values for residual
   * @param t Current time
   * @param fields Solution fields
   */
  void integrateResidual(const topology::Field& residual,
			 const PylithScalar t,
			 topology::SolutionFields* const fields);

  /** Integrate contributions to Jacobian matrix (A) associated with
   * operator.
   *
   * @param jacobian Sparse matrix for Jacobian of system.
   * @param t Current time
   * @param fields Solution fields
   */
  void integrateJacobian(topology::Jacobian* jacobian,
			 const PylithScalar t,
			 topology::SolutionFields* const fields);
  
// NOT IMPLEMENTED //////////////////////////////////////////////////////
private :

  /// Not implemented.
  ElasticityImplicitLgDeform(const ElasticityImplicitLgDeform&);

  /// Not implemented
  const ElasticityImplicitLgDeform& operator=(const ElasticityImplicitLgDeform&);

// PRIVATE MEMBERS //////////////////////////////////////////////////////
private :

  PylithScalar _dtm1; ///< Time step for t-dt1 -> t

}; // ElasticityImplicitLgDeform

#endif // pylith_feassemble_elasticityimplicitlgdeform_hh


// End of file 
