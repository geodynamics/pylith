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
 * @file pylith/feassemble/ExplicitElasticity.hh
 *
 * @brief Explicit time integration of dynamic elasticity equation
 * using finite-elements.
 *
 * Note: This object operates on a single finite-element family, which
 * is defined by the quadrature and a database of material property
 * parameters.
 *
 * Computes contributions to terms A and b in A(t) u(t+dt) = b(u(t),
 * u(t-dt)), where A(t) is a sparse matrix or vector, u(t+dt) is the
 * field we want to compute at time t+dt and b is a vector that
 * depends on the field at time t and t-dt.
 *
 * Contributions from elasticity include the intertial and stiffness
 * terms, so this object computes the following portions of A and b:
 *
 * A = 1/(dt*dt) [M]
 *
 * b = 2/(dt*dt)[M]{u(t)} - 1/(dt*dt)[M]{u(t-dt)} - [K]{u(t)}
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
 *   - Integrate and lump to form lumped matrix (field)
 *
 * See governing equations section of user manual for more
 * information.
 */

#if !defined(pylith_feassemble_explicitelasticity_hh)
#define pylith_feassemble_explicitelasticity_hh

#include "IntegratorExplicit.hh" // ISA IntegratorExplicit

namespace pylith {
  namespace feassemble {
    class ExplicitElasticity;
    class TestExplicitElasticity;
  } // feassemble

  namespace materials {
    class ElasticMaterial;
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

class pylith::feassemble::ExplicitElasticity : public IntegratorExplicit
{ // ExplicitElasticity
  friend class TestExplicitElasticity; // unit testing

// PUBLIC MEMBERS ///////////////////////////////////////////////////////
public :

  /// Constructor
  ExplicitElasticity(void);

  /// Destructor
  ~ExplicitElasticity(void);

  /// Create a copy of this object.
  IntegratorExplicit* clone(void) const;

  /** Set material.
   *
   * @param m Elastic material.
   */
  void material(const materials::ElasticMaterial* m);

  /** Integrate constant term (b) for dynamic elasticity term 
   * for 3-D finite elements.
   *
   * Compute b = 2/(dt*dt)[M]{u(t) - 1/(dt*dt)[M]{u(t-dt)} - [K]{u(t)}, where
   * [M] = density * [N]^T [N]
   *
   *
   * @param b Constant field (output)
   * @param dispT Displacement field at time t
   * @param dispTmdt Displacement field at time t-dt
   * @param coordinates Field of cell vertex coordinates
   */
  void integrateConstant(const ALE::Obj<real_section_type>& b,
			 const ALE::Obj<real_section_type>& dispT,
			 const ALE::Obj<real_section_type>& dispTmdt,
			 const ALE::Obj<real_section_type>& coordinates);

  /** Compute Jacobian matrix (A) associated with operator.
   *
   * @param mat Sparse matrix
   * @param dispT Displacement at time t
   * @param coordinates Field of cell vertex coordinates
   */
  void integrateJacobian(PetscMat* mat,
			 const ALE::Obj<real_section_type>& dispT,
			 const ALE::Obj<real_section_type>& coordinates);
  
// PROTECTED METHODS ////////////////////////////////////////////////////
protected :

  /** Copy constructor.
   *
   * @param i Integrator to copy
   */
  ExplicitElasticity(const ExplicitElasticity& i);

// PRIVATE METHODS //////////////////////////////////////////////////////
private :

  /// Not implemented
  const ExplicitElasticity& operator=(const ExplicitElasticity&);

// PRIVATE MEMBERS //////////////////////////////////////////////////////
private :

  /// Elastic material associated with integrator
  materials::ElasticMaterial* _material;

}; // ExplicitElasticity

#include "ExplicitElasticity.icc" // inline methods

#endif // pylith_feassemble_explicitelasticity_hh


// End of file 
