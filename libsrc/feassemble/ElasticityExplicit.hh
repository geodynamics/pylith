// -*- C++ -*-
//
// ======================================================================
//
//                           Brad T. Aagaard
//                        U.S. Geological Survey
//
// {LicenseText}
//
// ======================================================================
//

/**
 * @file pylith/feassemble/ElasticityExplicit.hh
 *
 * @brief Explicit time integration of dynamic elasticity equation
 * using finite-elements.
 *
 * Note: This object operates on a single finite-element family, which
 * is defined by the quadrature and a database of material property
 * parameters.
 *
 * Computes contributions to terms A and r in
 *
 * A(t) u(t+dt) = b(u(t), u(t-dt)),
 *
 * r = b - A u0(t+dt)
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

#if !defined(pylith_feassemble_elasticityexplicit_hh)
#define pylith_feassemble_elasticityexplicit_hh

#include "IntegratorElasticity.hh" // ISA IntegratorElasticity
#include "pylith/utils/array.hh" // USES std::vector, double_array

namespace pylith {
  namespace feassemble {
    class ElasticityExplicit;
  } // feassemble
} // pylith

namespace spatialdata {
  namespace spatialdb {
    class SpatialDB; // USES SpatialDB
  } // spatialdb
  namespace geocoords {
    class CoordSys; // USES CoordSys
  } // geocoords
} // spatialdata

class pylith::feassemble::ElasticityExplicit : public IntegratorElasticity
{ // ElasticityExplicit
  friend class TestElasticityExplicit; // unit testing

// PUBLIC MEMBERS ///////////////////////////////////////////////////////
public :

  /// Constructor
  ElasticityExplicit(void);

  /// Destructor
  ~ElasticityExplicit(void);

  /** Set time step for advancing from time t to time t+dt.
   *
   * @param dt Time step
   */
  void timeStep(const double dt);

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
   * @param mesh Finite-element mesh
   */
  void integrateResidual(const ALE::Obj<real_section_type>& residual,
			 const double t,
			 topology::FieldsManager* const fields,
			 const ALE::Obj<Mesh>& mesh,
			 const spatialdata::geocoords::CoordSys* cs);

  /** Integrate contributions to Jacobian matrix (A) associated with
   * operator.
   *
   * @param jacobian Sparse matrix to hold Jacobian of operator.
   * @param t Current time
   * @param fields Solution fields.
   * @param mesh Finite-element mesh.
   */
  void integrateJacobian(PetscMat* jacobian,
			 const double t,
			 topology::FieldsManager* const fields,
			 const ALE::Obj<Mesh>& mesh);

// NOT IMPLEMENTED //////////////////////////////////////////////////////
private :

  /// Not implemented.
  ElasticityExplicit(const ElasticityExplicit& i);

  /// Not implemented
  const ElasticityExplicit& operator=(const ElasticityExplicit&);

// PRIVATE MEMBERS //////////////////////////////////////////////////////
private :

  double _dtm1; ///< Time step for t-dt1 -> t

}; // ElasticityExplicit

#endif // pylith_feassemble_elasticityexplicit_hh


// End of file 
