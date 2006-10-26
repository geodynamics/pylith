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
 * @file pylith/feassemble/IntegratorInertia3D.hh
 *
 * @brief Integrate inertial term for 3-D finite elements.
 */

#if !defined(pylith_feassemble_integratorinertia3d_hh)
#define pylith_feassemble_integratorinertia3d_hh

#include "Integrator.hh"

namespace pylith {
  namespace feassemble {
    class IntegratorInertia3D;
    class TestIntegratorInertia3D;
  } // feassemble
} // pylith

class pylith::feassemble::IntegratorInertia3D : public Integrator
{ // Integrator1D
  friend class TestIntegratorInertia3D; // unit testing

// PUBLIC MEMBERS ///////////////////////////////////////////////////////
public :

  /// Constructor
  IntegratorInertia3D(void);

  /// Destructor
  ~IntegratorInertia3D(void);

  /// Create a copy of this object.
  Integrator* clone(void) const;

  /** Integrate inertial term for 3-D finite elements.
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
		 const ALE::Obj<ALE::Mesh::real_section_type>& coordinates);

// PROTECTED METHODS ////////////////////////////////////////////////////
protected :

  /** Copy constructor.
   *
   * @param i Integrator to copy
   */
  IntegratorInertia3D(const IntegratorInertia3D& i);

// PRIVATE METHODS //////////////////////////////////////////////////////
private :

  /// Not implemented
  const IntegratorInertia3D& operator=(const IntegratorInertia3D&);

}; // IntegratorInertia3D

#include "IntegratorInertia3D.icc" // inline methods

#endif // pylith_feassemble_integratorineria3d_hh

// End of file 
