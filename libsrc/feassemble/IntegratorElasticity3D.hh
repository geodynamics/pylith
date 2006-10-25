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
   * @param 
   */
  void integrateAction(const ALE::Mesh::Obj<section_type>& A,
		       const ALE::Mesh::Obj<section_type>& F,
		       const ALE::Mesh::Obj<section_type>& coordinates);

  /** Compute tangent stiffness matrix.
   *
   * @param
   */
  void integrate(PetscMat& mat,
		 const ALE::Mesh::Obj<section_type>& A,
		 const ALE::Mesh::Obj<section_type>& F,
		 const ALE::Mesh::Obj<section_type>& coordinates);

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

#endif // pylith_feassemble_quadrature3d_hh

// End of file 
