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
 * @file libsrc/feassemble/ElasticityExplicitLgDeform.hh
 *
 * @brief Explicit time integration of dynamic elasticity equation
 * with large rigid body motion and small strains using
 * finite-elements.
 */

#if !defined(pylith_feassemble_elasticityexplicitlgdeform_hh)
#define pylith_feassemble_elasticityexplicitlgdeform_hh

// Include directives ---------------------------------------------------
#include "IntegratorElasticityLgDeform.hh" // ISA IntegratorElasticityLgDeform

// ElasticityExplicitLgDeform -------------------------------------------
/**@brief Explicit time integration of the dynamic elasticity equation
 * with large rigid body motion and small strains using finite-elements.
 *
 * Note: This object operates on a single finite-element family, which
 * is defined by the quadrature and a database of material property
 * parameters.
 *
 * Computes contributions to terms A and r in
 *
 * A(t+dt) du(t) = b(t+dt, u(t), u(t-dt)) - A(t+dt) u(t),
 *
 * r(t+dt) = b(t+dt) - A(t+dt) (u(t) + du(t))
 *
 * where A(t) is a sparse matrix or vector, u(t+dt) is the field we
 * want to compute at time t+dt, b is a vector that depends on the
 * field at time t and t-dt, and u0 is zero at unknown DOF and set to
 * the known values at the constrained DOF.
 *
 * Contributions from elasticity include the intertial and stiffness
 * terms, so this object computes the following portions of A and r:
 *
 * A = 1/(dt*dt) [M]
 *
 * r = (1/(dt*dt) [M])(- {u(t+dt)} + 2/(dt*dt){u(t)} - {u(t-dt)}) - [K]{u(t)}
 *
 * Translational inertia.
 *   - Residual action over cell
 *     \f[
 *       \int_{V^e} \rho N^p \sum_q N^q u_i^q \, dV
 *     \f]
 *   - Jacobian action over cell
 *     \f[
 *       \int_{V^e} (\rho N^q N^q)_i \, dV
 *     \f]
 *
 * See governing equations section of user manual for more
 * information.
*/
class pylith::feassemble::ElasticityExplicitLgDeform : public IntegratorElasticityLgDeform
{ // ElasticityExplicitLgDeform
  friend class TestElasticityExplicitLgDeform; // unit testing

// PUBLIC MEMBERS ///////////////////////////////////////////////////////
public :

  /// Constructor
  ElasticityExplicitLgDeform(void);

  /// Destructor
  ~ElasticityExplicitLgDeform(void);

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

  /** Set normalized viscosity for numerical damping.
   *
   * @param viscosity Normalized viscosity (viscosity / elastic modulus).
   */
  void normViscosity(const PylithScalar viscosity);

  /** Integrate contributions to residual term (r) for operator.
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
   * @param jacobian Diagonal matrix (as field) for Jacobian of system.
   * @param t Current time
   * @param fields Solution fields
   */
  void integrateJacobian(topology::Field* jacobian,
			 const PylithScalar t,
			 topology::SolutionFields* const fields);

// NOT IMPLEMENTED //////////////////////////////////////////////////////
private :

  /// Not implemented.
  ElasticityExplicitLgDeform(const ElasticityExplicitLgDeform&);

  /// Not implemented
  const ElasticityExplicitLgDeform& operator=(const ElasticityExplicitLgDeform&);

  /// Not implemented.
  void integrateJacobian(topology::Jacobian*,
			 const PylithScalar,
			 topology::SolutionFields* const);


// PRIVATE MEMBERS //////////////////////////////////////////////////////
private :

  PylithScalar _dtm1; ///< Time step for t-dt1 -> t
  PylithScalar _normViscosity; ///< Normalized viscosity for numerical damping.

}; // ElasticityExplicitLgDeform

#endif // pylith_feassemble_elasticityexplicitlgdeform_hh


// End of file 
