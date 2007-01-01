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
 * @file pylith/feassemble/IntegratorElasticity3D.hh
 *
 * @brief Integrate elasticity term for 3-D finite elements.
 */

#if !defined(pylith_feassemble_integratorelasticity3d_hh)
#define pylith_feassemble_integratorelasticity3d_hh

#include "Integrator.hh"

namespace pylith {
  namespace feassemble {
    class IntegratorElasticity3D;
    class TestIntegratorElasticity3D;
  } // feassemble
} // pylith

class pylith::feassemble::IntegratorElasticity3D : public Integrator
{ // Integrator1D
  friend class TestIntegratorElasticity3D; // unit testing

// PUBLIC MEMBERS ///////////////////////////////////////////////////////
public :

  /// Constructor
  IntegratorElasticity3D(void);

  /// Destructor
  ~IntegratorElasticity3D(void);

  /// Create a copy of this object.
  Integrator* clone(void) const;

  /** Integrate elasticity term for 3-D finite elements.
   *
   * @param fieldOut Output field
   * @param fieldIn Input field
   * @param coordinates Field of cell vertex coordinates
   */
  void integrateAction(const ALE::Obj<ALE::Mesh::real_section_type>& fieldOut,
		       const ALE::Obj<ALE::Mesh::real_section_type>& fieldIn,
		       const ALE::Obj<ALE::Mesh::real_section_type>& coordinates);

  /** Compute matrix associated with operator.
   *
   * @param mat Sparse matrix
   * @param coordinates Field of cell vertex coordinates
   */
  void integrate(PetscMat* mat,
		 const ALE::Obj<ALE::Mesh::real_section_type>& fieldIn,
		 const ALE::Obj<ALE::Mesh::real_section_type>& coordinates);

  /** Initialize, get material property parameters from database.
   *
   * @param mesh PETSc mesh
   * @param db Pointer to spatial database with material property parameters
   * @param cs Pointer to coordinate system of vertices
   */
  void initialize(ALE::Obj<ALE::Mesh>& mesh,
		  spatialdata::spatialdb::SpatialDB* db,
		  spatialdata::geocoords::CoordSys* cs);

// PROTECTED METHODS ////////////////////////////////////////////////////
protected :

  /** Copy constructor.
   *
   * @param i Integrator to copy
   */
  IntegratorElasticity3D(const IntegratorElasticity3D& i);

// PRIVATE METHODS //////////////////////////////////////////////////////
private :

  /// Not implemented
  const IntegratorElasticity3D& operator=(const IntegratorElasticity3D&);

}; // IntegratorElasticity3D

#include "IntegratorElasticity3D.icc" // inline methods

#endif // pylith_feassemble_integratorelasticity3d_hh

// End of file 
