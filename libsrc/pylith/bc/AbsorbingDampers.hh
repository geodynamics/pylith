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
// Copyright (c) 2010-2015 University of California, Davis
//
// See COPYING for license information.
//
// ----------------------------------------------------------------------
//

/** @file libsrc/bc/AbsorbingDampers.hh
 *
 * @brief C++ implementation of an absorbing boundary with simple
 * dampers.
 *
 * Computes contributions to terms A and r in 
 *
 * A(t) u(t+dt) = b(u(t), u(t-dt)),
 *
 * r = b - a u0(t+dt)
 *
 * where A(t) is a sparse matrix or vector, u(t+dt) is the field we
 * want to compute at time t+dt, b is a vector that depends on the
 * field at time t and t-dt, and u0 is zero at unknown DOF and set to
 * the known values at the constrained DOF.
 *
 * Contributions from elasticity include the intertial and stiffness
 * terms, so this object computes the following portions of A and r:
 *
 * A = 1/(2*dt) [C]
 *
 * r = (1/(2*dt) [C])(- {u(t+dt)} + {u(t-dt)})
 *
 * Damping matrix.
 *   - Residual action over cell
 *     \f[
 *       \int_{S^e} \rho c_i N^p \sum_q N^q u_i^q \, dS
 *     \f]
 *   - Jacobian action over cell
 *     \f[
 *       \int_{S^e} (\rho c_i N^q N^q)_i \, dS
 *     \f]
 *
 * See boundary conditions section of user manual for more
 * information.
 */

#if !defined(pylith_bc_absorbingdampers_hh)
#define pylith_bc_absorbingdampers_hh

// Include directives ---------------------------------------------------
#include "BCIntegratorSubMesh.hh" // ISA BCIntegratorSubMesh

// AbsorbingDampers ------------------------------------------------------
/// Absorbing boundary with simple dampers.
class pylith::bc::AbsorbingDampers : public BCIntegratorSubMesh
{ // class AbsorbingDampers
  friend class TestAbsorbingDampers; // unit testing

  // PUBLIC METHODS /////////////////////////////////////////////////////
public :

  /// Default constructor.
  AbsorbingDampers(void);

  /// Destructor.
  ~AbsorbingDampers(void);

  /// Deallocate PETSc and local data structures.
  void deallocate(void);
  
  /** Set database for boundary condition parameters.
   *
   * @param db Spatial database
   */
  void db(spatialdata::spatialdb::SpatialDB* const db);

  /** Initialize boundary condition.
   *
   * @param mesh Finite-element mesh.
   * @param upDir Direction perpendicular to horizontal surface tangent 
   *   direction that is not collinear with surface normal.
   */
  void initialize(const topology::Mesh& mesh,
		  const PylithScalar upDir[3]);

  /** Integrate contributions to residual term (r) for operator.
   *
   * @param residual Field containing values for residual
   * @param t Current time
   * @param fields Solution fields
   */
  void integrateResidual(const topology::Field& residual,
			 const PylithScalar t,
			 topology::SolutionFields* const fields);

  /** Integrate contributions to residual term (r) for operator.
   *
   * @param residual Field containing values for residual
   * @param t Current time
   * @param fields Solution fields
   */
  void integrateResidualLumped(const topology::Field& residual,
			       const PylithScalar t,
			       topology::SolutionFields* const fields);

  /** Integrate contributions to Jacobian matrix (A) associated with
   * operator.
   *
   * @param jacobian Jacobian of system.
   * @param t Current time
   * @param fields Solution fields
   */
  void integrateJacobian(topology::Jacobian* jacobian,
			 const PylithScalar t,
			 topology::SolutionFields* const fields);

  /** Integrate contributions to Jacobian matrix (A) associated with
   * operator.
   *
   * @param jacobian Jacobian of system.
   * @param t Current time
   * @param fields Solution fields
   */
  void integrateJacobian(topology::Field* jacobian,
			 const PylithScalar t,
			 topology::SolutionFields* const fields);

  /** Verify configuration is acceptable.
   *
   * @param mesh Finite-element mesh
   */
  void verifyConfiguration(const topology::Mesh& mesh) const;

  // PRIVATE METHODS ////////////////////////////////////////////////////
private :

  /// Initialize logger.
  void _initializeLogger(void);

  // PRIVATE MEMBERS ////////////////////////////////////////////////////
private :

  topology::VecVisitorSubMesh* _velocityVisitor; ///< Cache velocity field visitor.

  spatialdata::spatialdb::SpatialDB* _db; ///< Spatial database w/parameters

  // NOT IMPLEMENTED ////////////////////////////////////////////////////
private :

  /// Not implemented
  AbsorbingDampers(const AbsorbingDampers&);

  /// Not implemented
  const AbsorbingDampers& operator=(const AbsorbingDampers&);

}; // class AbsorbingDampers

#include "AbsorbingDampers.icc" // inline methods

#endif // pylith_bc_absorbingdampers_hh


// End of file 
