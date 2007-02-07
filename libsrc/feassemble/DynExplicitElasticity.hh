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
 * @file pylith/feassemble/DynExplicitElasticity.hh
 *
 * @brief Explicit time integration of dynamic elasticity equation
 * using finite-elements.
 *
 * Note: This object operates on a single finite-element family, which
 * is defined by the quadrature and a database of material property
 * parameters.
 *
 * Computes terms A and b in A(t) u(t+dt) = b(u(t), u(t-dt)), where
 * A(t) is a sparse matrix or vector, u(t+dt) is the field we want to
 * compute at time t+dt and b is a vector that depends on the field at
 * time t and t-dt.
 *
 * A = 1/(dt*dt) [M]
 *
 * b = 2/(dt*dt)[M]{u(t)} - 1/(dt*dt)[M]{u(t-dt)} - [K]{u(t)} + {f(t)}
 */

#if !defined(pylith_feassemble_dynexplicitelasticity_hh)
#define pylith_feassemble_dynexplicitelasticity_hh

#include "IntegratorDynExplicit.hh" // ISA IntegratorDynExplicit

namespace pylith {
  namespace feassemble {
    class DynExplicitElasticity;
    class TestDynExplicitElasticity;
  } // feassemble
} // pylith

class pylith::feassemble::DynExplicitElasticity : public IntegratorDynExplicit
{ // DynExplicitElasticity
  friend class TestDynExplicitElasticity; // unit testing

// PUBLIC MEMBERS ///////////////////////////////////////////////////////
public :

  /// Constructor
  DynExplicitElasticity(void);

  /// Destructor
  ~DynExplicitElasticity(void);

  /// Create a copy of this object.
  IntegratorDynExplicit* clone(void) const;

  /** Integrate residual term (b) for dynamic elasticity term 
   * for 3-D finite elements.
   *
   * @param residual Residual field (output)
   * @param dispT Displacement field at time t
   * @param dispTmdt Displacement field at time t-dt
   * @param coordinates Field of cell vertex coordinates
   */
  void integrateResidual(const ALE::Obj<real_section_type>& residual,
			 const ALE::Obj<real_section_type>& dispT,
			 const ALE::Obj<real_section_type>& dispTmdt,
			 const ALE::Obj<real_section_type>& coordinates);

  /** Compute matrix (A) associated with operator.
   *
   * @param mat Sparse matrix
   * @param dispT Displacement at time t
   * @param coordinates Field of cell vertex coordinates
   */
  void integrateJacobian(PetscMat* mat,
			 const ALE::Obj<real_section_type>& dispT,
			 const ALE::Obj<real_section_type>& coordinates);
  
  /** Compute field (A) associated with lumped operator.
   *
   * @param fieldOut Output Jacobian field
   * @param dispT Displacement at time t
   * @param coordinates Field of cell vertex coordinates
   */
  void integrateJacobian(const ALE::Obj<real_section_type>& fieldOut,
			 const ALE::Obj<real_section_type>& dispT,
			 const ALE::Obj<real_section_type>& coordinates);
  
  /** Initialize, get material property parameters from database.
   *
   * @param mesh PETSc mesh
   * @param cs Pointer to coordinate system of vertices
   * @param db Pointer to spatial database with material property parameters
   */
  void initialize(ALE::Obj<ALE::Mesh>& mesh,
		  spatialdata::geocoords::CoordSys* cs,
		  spatialdata::spatialdb::SpatialDB* db);

// PROTECTED METHODS ////////////////////////////////////////////////////
protected :

  /** Copy constructor.
   *
   * @param i Integrator to copy
   */
  DynExplicitElasticity(const DynExplicitElasticity& i);

// PRIVATE METHODS //////////////////////////////////////////////////////
private :

  /// Not implemented
  const DynExplicitElasticity& operator=(const DynExplicitElasticity&);

}; // DynExplicitElasticity

#include "DynExplicitElasticity.icc" // inline methods

#endif // pylith_feassemble_dynexplicitelasticity_hh


// End of file 
