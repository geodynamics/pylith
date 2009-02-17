// -*- C++ -*-
//
// ----------------------------------------------------------------------
//
//                           Brad T. Aagaard
//                        U.S. Geological Survey
//
// {LicenseText}
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
#include "BoundaryCondition.hh" // ISA BoundaryCondition

#include "pylith/topology/SubMesh.hh" // ISA Quadrature<SubMesh>
#include "pylith/feassemble/Quadrature.hh" // ISA Integrator<Quadrature>
#include "pylith/feassemble/Integrator.hh" // ISA Integrator

#include "pylith/utils/array.hh" // USES std::vector, double_array, int_array

// AbsorbingDampers ------------------------------------------------------
class pylith::bc::AbsorbingDampers : public BoundaryCondition, 
				     public feassemble::Integrator<feassemble::Quadrature<topology::SubMesh> >
{ // class AbsorbingDampers
  friend class TestAbsorbingDampers; // unit testing

  // PUBLIC METHODS /////////////////////////////////////////////////////
public :

  /// Default constructor.
  AbsorbingDampers(void);

  /// Destructor.
  ~AbsorbingDampers(void);

  /** Initialize boundary condition.
   *
   * @param mesh Finite-element mesh.
   * @param upDir Direction perpendicular to horizontal surface tangent 
   *   direction that is not collinear with surface normal.
   */
  void initialize(const topology::Mesh& mesh,
		  const double upDir[3]);

  /** Integrate contributions to residual term (r) for operator.
   *
   * @param residual Field containing values for residual
   * @param t Current time
   * @param fields Solution fields
   */
  void integrateResidual(const topology::Field<topology::Mesh>& residual,
			 const double t,
			 topology::SolutionFields* const fields);

  /** Integrate contributions to Jacobian matrix (A) associated with
   * operator.
   *
   * @param mat Sparse matrix
   * @param t Current time
   * @param fields Solution fields
   */
  void integrateJacobian(PetscMat* mat,
			 const double t,
			 topology::SolutionFields* const fields);

  /** Verify configuration is acceptable.
   *
   * @param mesh Finite-element mesh
   */
  void verifyConfiguration(const topology::Mesh& mesh) const;

  // NOT IMPLEMENTED ////////////////////////////////////////////////////
private :

  /// Not implemented
  AbsorbingDampers(const AbsorbingDampers& m);

  /// Not implemented
  const AbsorbingDampers& operator=(const AbsorbingDampers& m);

  // PRIVATE MEMBERS ////////////////////////////////////////////////////
private :

  /// Mesh of absorbing boundary
  topology::SubMesh* _boundaryMesh;

  /// Damping constants in global coordinates at integration points.
  topology::Field<topology::SubMesh>* _dampingConsts;

}; // class AbsorbingDampers

#endif // pylith_bc_absorbingdampers_hh


// End of file 
