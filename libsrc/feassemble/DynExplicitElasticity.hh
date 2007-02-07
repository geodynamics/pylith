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
 * Computes terms A and b in A(t) u(t+dt) = b(u(t), u(t-dt)), where
 * A(t) is a sparse matrix or vector, u(t+dt) is the field we want to
 * compute at time t+dt and b is a vector that depends on the field at
 * time t and t-dt.
 *
 * A = 1/(dt*dt) [M]
 *
 * b = 2/(dt*dt)[M]{u(t)} - 1/(dt*dt)[M]{u(t-dt)} - 1/s[K]{u(t)} + {f(t)}
 */

#if !defined(pylith_feassemble_dynexplicitelasticity_hh)
#define pylith_feassemble_dynexplicitelasticity_hh

#include <petscmesh.h> // USES Mesh
#include "pylith/utils/petscfwd.h" // USES PetscMat

#include "IntegratorDynExplicit.hh"

namespace pylith {
  namespace feassemble {
    class DynExplicitElasticity;
    class TestDynExplicitElasticity;

    class Quadrature; // HOLDSA Quadrature
  } // feassemble
} // pylith

class pylith::feassemble::DynExplicitElasticity : 
  public IntegratorDynExplicit
{ // DynExplicitElasticity
  friend class TestDynExplicitElasticity; // unit testing

// PUBLIC TYPEDEFS //////////////////////////////////////////////////////
public :

  typedef ALE::Mesh Mesh;
  typedef Mesh::topology_type topology_type;
  typedef topology_type::point_type point_type;
  typedef Mesh::real_section_type real_section_type;

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
   * @param fieldOut Output field
   * @param fieldInT Input field at time t
   * @param fieldInTmdt Input field at time t-dt
   * @param coordinates Field of cell vertex coordinates
   */
  void integrateResidual(const ALE::Obj<real_section_type>& fieldOut,
			 const ALE::Obj<real_section_type>& fieldInT,
			 const ALE::Obj<real_section_type>& fieldInTmdt,
			 const ALE::Obj<real_section_type>& coordinates);

  /** Compute matrix (A) associated with operator.
   *
   * @param mat Sparse matrix
   * @param fieldIn Input field at time t
   * @param coordinates Field of cell vertex coordinates
   */
  void integrateJacobian(PetscMat* mat,
			 const ALE::Obj<real_section_type>& fieldIn,
			 const ALE::Obj<real_section_type>& coordinates);
  
  /** Compute field (A) associated with lumped operator.
   *
   * @param fieldOut Output Jacobian field
   * @param fieldIn Input field at time t
   * @param coordinates Field of cell vertex coordinates
   */
  void integrateJacobian(const ALE::Obj<real_section_type>& fieldOut,
			 const ALE::Obj<real_section_type>& fieldIn,
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
