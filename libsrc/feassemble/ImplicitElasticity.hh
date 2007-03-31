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
 * @file pylith/feassemble/ImplicitElasticity.hh
 *
 * @brief Implicit time integration of quasi-static elasticity equation
 * using finite-elements.
 *
 * Note: This object operates on a single finite-element family, which
 * is defined by the quadrature and a database of material property
 * parameters.
 *
 * Computes contributions to terms A and b in A(t+dt) u(t+dt) = b(u(t),
 * u(t+dt)), where A(t+dt) is a sparse matrix, u(t+dt) is the
 * field we want to compute at time t+dt and b is a vector that
 * depends on the field at time t and t+dt.
 *
 * Contributions to the RHS (b) include body forces, which are either
 * independent of u (small strain case) or are computed based on the
 * displacements at time t. The RHS also includes the internal force
 * vector, which is either constant for a time step (small strain,
 * linear rheology) or changes with each iteration (large strain or
 * nonlinear rheology). The internal force vector is subtracted from the
 * existing force vector to get the residual. This object also computes
 * the entire stiffness matrix (A).
 *
 * See governing equations section of user manual for more
 * information.
 */

#if !defined(pylith_feassemble_implicitelasticity_hh)
#define pylith_feassemble_implicitelasticity_hh

#include "IntegratorImplicit.hh" // ISA IntegratorImplicit

namespace pylith {
  namespace feassemble {
    class ImplicitElasticity;
    class TestImplicitElasticity;
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

class pylith::feassemble::ImplicitElasticity : public IntegratorImplicit
{ // ImplicitElasticity
  friend class TestImplicitElasticity; // unit testing

// PUBLIC MEMBERS ///////////////////////////////////////////////////////
public :

  /// Constructor
  ImplicitElasticity(void);

  /// Destructor
  ~ImplicitElasticity(void);

  /// Create a copy of this object.
  IntegratorImplicit* clone(void) const;

  /** Set material.
   *
   * @param m Elastic material.
   */
  void material(const materials::ElasticMaterial* m);

  /** Integrate constant part of RHS for 3-D finite elements.
   * Includes only body forces at present.
   *
   *
   * @param b Constant field (output)
   * @param dispT Displacement field at time t
   * @param grav Gravity vector for cell
   * @param coordinates Field of cell vertex coordinates
   */
  void integrateConstant(const ALE::Obj<real_section_type>& b,
			 const ALE::Obj<real_section_type>& dispT,
			 const ALE::Obj<real_section_type>& grav,
			 const ALE::Obj<real_section_type>& coordinates);

  /** Integrate residual part of RHS for 3-D finite elements.
   * Includes element internal force contribution.
   *
   *
   * @param b Constant field (output)
   * @param dispT Displacement field at time t
   * @param coordinates Field of cell vertex coordinates
   */
  void integrateConstant(const ALE::Obj<real_section_type>& b,
			 const ALE::Obj<real_section_type>& dispT,
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
  ImplicitElasticity(const ImplicitElasticity& i);

// PRIVATE METHODS //////////////////////////////////////////////////////
private :

  /// Not implemented
  const ImplicitElasticity& operator=(const ImplicitElasticity&);

// PRIVATE MEMBERS //////////////////////////////////////////////////////
private :

  /// Elastic material associated with integrator
  materials::ElasticMaterial* _material;

}; // ImplicitElasticity

#include "ImplicitElasticity.icc" // inline methods

#endif // pylith_feassemble_implicitelasticity_hh


// End of file 
