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
 * @file pylith/feassemble/Integrator.hh
 *
 * @brief Abstract base class for integration of finite-element actions.
 */

#if !defined(pylith_feassemble_integrator_hh)
#define pylith_feassemble_integrator_hh

#include <Mesh.hh> // USES Mesh
#include "pylith/utils/petscfwd.h" // USES PetscMat

namespace pylith {
  namespace feassemble {
    class Integrator;
    class TestIntegrator;

    class Quadrature; // HOLDSA Quadrature
  } // feassemble
} // pylith

class pylith::feassemble::Integrator
{ // Integrator
  friend class TestIntegrator; // unit testing

// PUBLIC MEMBERS ///////////////////////////////////////////////////////
public :

  /// Constructor
  Integrator(void);

  /// Destructor
  virtual
  ~Integrator(void);

  /// Create a copy of this object.
  virtual
  Integrator* clone(void) const = 0;

  /** Integrate elasticity term for 3-D finite elements.
   *
   * @param fieldOut Output field
   * @param fieldIn Input field
   * @param coordinates Field of cell vertex coordinates
   */
  virtual 
  void integrateAction(const ALE::Obj<ALE::Mesh::real_section_type>& fieldOut,
		       const ALE::Obj<ALE::Mesh::real_section_type>& fieldIn,
		       const ALE::Obj<ALE::Mesh::real_section_type>& coordinates) = 0;

  /** Compute matrix associated with operator.
   *
   * @param mat Sparse matrix
   * @param coordinates Field of cell vertex coordinates
   */
  virtual 
  void integrate(PetscMat* mat,
		 const ALE::Obj<ALE::Mesh::real_section_type>& coordinates) = 0;

  /** Set quadrature for integrating finite-element quantities.
   *
   * @param q Quadrature for integrating.
   */
  void quadrature(const Quadrature* q);

// PROTECTED METHODS ////////////////////////////////////////////////////
protected :

  /** Copy constructor.
   *
   * @param i Integrator to copy
   */
  Integrator(const Integrator& i);

// PRIVATE METHODS //////////////////////////////////////////////////////
private :

  /// Not implemented
  const Integrator& operator=(const Integrator&);

// PRIVATE MEMBERS //////////////////////////////////////////////////////
private :

  Quadrature* _quadrature; ///< Quadrature for integrating finite-element

}; // Integrator

#endif // pylith_feassemble_integrator_hh

// End of file 
