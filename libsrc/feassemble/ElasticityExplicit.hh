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
// Copyright (c) 2010 University of California, Davis
//
// See COPYING for license information.
//
// ======================================================================
//

/**
 * @file libsrc/feassemble/ElasticityExplicit.hh
 *
 * @brief Explicit time integration of dynamic elasticity equation
 * using finite-elements.
 */

#if !defined(pylith_feassemble_elasticityexplicit_hh)
#define pylith_feassemble_elasticityexplicit_hh

// Include directives ---------------------------------------------------
#include "IntegratorElasticity.hh" // ISA IntegratorElasticity

// ElasticityExplicit ---------------------------------------------------
/**@brief Explicit time integration of the dynamic elasticity equation
 * using finite-elements.
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
class pylith::feassemble::ElasticityExplicit : public IntegratorElasticity
{ // ElasticityExplicit
  friend class TestElasticityExplicit; // unit testing

// PUBLIC MEMBERS ///////////////////////////////////////////////////////
public :

  /// Constructor
  ElasticityExplicit(void);

  /// Destructor
  ~ElasticityExplicit(void);

  /// Deallocate PETSc and local data structures.
  void deallocate(void);
  
  /** Set time step for advancing from time t to time t+dt.
   *
   * @param dt Time step
   */
  void timeStep(const double dt);

  /** Set normalized viscosity for numerical damping.
   *
   * @param viscosity Normalized viscosity (viscosity / elastic modulus).
   */
  void normViscosity(const double viscosity);

  /** Set flag for setting constraints for total field solution or
   *  incremental field solution.
   *
   * @param flag True if using incremental solution, false otherwise.
   */
  void useSolnIncr(const bool flag);

  /** Integrate contributions to residual term (r) for operator.
   *
   * @param residual Field containing values for residual
   * @param t Current time
   * @param fields Solution fields
   */
  void integrateResidual(const topology::Field<topology::Mesh>& residual,
			 const double t,
			 topology::SolutionFields* const fields);

  /** Integrate contributions to residual term (r) for operator.
   *
   * @param residual Field containing values for residual
   * @param t Current time
   * @param fields Solution fields
   */
  void integrateResidualLumped(const topology::Field<topology::Mesh>& residual,
       const double t,
       topology::SolutionFields* const fields);

  /** Integrate contributions to Jacobian matrix (A) associated with
   * operator.
   *
   * @param jacobian Sparse matrix for Jacobian of system.
   * @param t Current time
   * @param fields Solution fields
   */
  void integrateJacobian(topology::Jacobian* jacobian,
			 const double t,
			 topology::SolutionFields* const fields);

  /** Integrate contributions to Jacobian matrix (A) associated with
   * operator.
   *
   * @param jacobian Diagonal matrix (as field) for Jacobian of system.
   * @param t Current time
   * @param fields Solution fields
   */
  void integrateJacobian(topology::Field<topology::Mesh>* jacobian,
			 const double t,
			 topology::SolutionFields* const fields);

// NOT IMPLEMENTED //////////////////////////////////////////////////////
private :

  /// Not implemented.
  ElasticityExplicit(const ElasticityExplicit&);

  /// Not implemented
  const ElasticityExplicit& operator=(const ElasticityExplicit&);

// PRIVATE MEMBERS //////////////////////////////////////////////////////
private :

  double _dtm1; ///< Time step for t-dt1 -> t
  double _normViscosity; ///< Normalized viscosity for numerical damping.

}; // ElasticityExplicit

#endif // pylith_feassemble_elasticityexplicit_hh


// End of file 
