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
 * @file pylith/feassemble/ElasticityImplicit.hh
 *
 * @brief Implicit time integration of quasi-static elasticity equation
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

#if !defined(pylith_feassemble_elasticityimplicit_hh)
#define pylith_feassemble_elasticityimplicit_hh

#include "Integrator.hh" // ISA Integrator

namespace pylith {
  namespace feassemble {
    class ElasticityImplicit;
    class TestElasticityImplicit;
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

class pylith::feassemble::ElasticityImplicit : public Integrator
{ // ElasticityImplicit
  friend class TestElasticityImplicit; // unit testing

// PUBLIC MEMBERS ///////////////////////////////////////////////////////
public :

  /// Constructor
  ElasticityImplicit(void);

  /// Destructor
  ~ElasticityImplicit(void);

  /** Set material.
   *
   * @param m Elastic material.
   */
  void material(materials::ElasticMaterial* m);

  /** Set time step for advancing from time t to time t+dt.
   *
   * @param dt Time step
   */
  void timeStep(const double dt);

  /** Integrate residual part of RHS for 3-D finite elements.
   * Includes gravity and element internal force contribution.
   *
   * We assume that the effects of boundary conditions are already
   * included in b (tractions, concentrated nodal forces, and
   * contributions to internal force vector due to
   * displacement/velocity BC).  This routine computes the additional
   * external loads due to body forces (not yet implemented) plus the
   * element internal forces for the current stress state.
   *
   * @param b Residual field (output)
   * @param fields Solution fields
   * @param mesh Mesh object
   */
  void integrateResidual(const ALE::Obj<real_section_type>& b,
			 topology::FieldsManager* const fields,
			 const ALE::Obj<Mesh>& mesh);

  /** Compute Jacobian matrix (A) associated with operator.
   *
   * @param mat Sparse matrix
   * @param fields Solution fields
   * @param mesh Mesh object
   */
  void integrateJacobian(PetscMat* mat,
			 topology::FieldsManager* const fields,
			 const ALE::Obj<Mesh>& mesh);
  
  /** Update state variables as needed.
   *
   * @param field Current solution field.
   * @param mesh Finite-element mesh
   */
  void updateState(const ALE::Obj<real_section_type>& field,
		   const ALE::Obj<Mesh>& mesh);

// PRIVATE METHODS //////////////////////////////////////////////////////
private :

  /// Not implemented.
  ElasticityImplicit(const ElasticityImplicit& i);

  /// Not implemented
  const ElasticityImplicit& operator=(const ElasticityImplicit&);

// PRIVATE MEMBERS //////////////////////////////////////////////////////
private :

  double _dtm1; ///< Time step for t-dt1 -> t

  /// Elastic material associated with integrator
  materials::ElasticMaterial* _material;

}; // ElasticityImplicit

#endif // pylith_feassemble_elasticityimplicit_hh


// End of file 
