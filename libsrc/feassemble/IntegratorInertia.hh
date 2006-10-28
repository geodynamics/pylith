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
 * @file pylith/feassemble/IntegratorInertia.hh
 *
 * @brief Integrate inertial term for 3-D finite elements.
 *
 * This object performs integrals over the domain of a finite-element
 * associated with translational inertia.
 *   - Integrate action over cell
 *     \f[
 *       \int_{V^e} \rho N^p \sum_q N^q u_i^q \, dV
 *     \f]
 *   - Integrate to form matrix
 *     \f[
 *       \int_{V^e} (\rho N^q N^q)_i \, dV
 *     \f]
 *   - Integrate and lump to form lumped matrix (field)
 *
 * See governing equations section of user manual for more
 * information.
 */

#if !defined(pylith_feassemble_integratorinertia_hh)
#define pylith_feassemble_integratorinertia_hh

#include "Integrator.hh"

namespace pylith {
  namespace feassemble {
    class IntegratorInertia;
    class TestIntegratorInertia;
  } // feassemble
} // pylith

class pylith::feassemble::IntegratorInertia : public Integrator
{ // Integrator1D
  friend class TestIntegratorInertia; // unit testing

// PUBLIC MEMBERS ///////////////////////////////////////////////////////
public :

  /// Constructor
  IntegratorInertia(void);

  /// Destructor
  ~IntegratorInertia(void);

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
		 const ALE::Obj<real_section_type>& coordinates);

  /** Compute lumped mass matrix.
   *
   * @param fieldOut Lumped mass matrix as field
   * @param coordinates Field of cell vertex coordinates
   */
  void integrateLumped(const ALE::Obj<real_section_type>& fieldOut,
		       const ALE::Obj<real_section_type>& coordinates);

// PROTECTED METHODS ////////////////////////////////////////////////////
protected :

  /** Copy constructor.
   *
   * @param i Integrator to copy
   */
  IntegratorInertia(const IntegratorInertia& i);

// PRIVATE METHODS //////////////////////////////////////////////////////
private :

  /// Not implemented
  const IntegratorInertia& operator=(const IntegratorInertia&);

}; // IntegratorInertia

#include "IntegratorInertia.icc" // inline methods

#endif // pylith_feassemble_integratorineria_hh

// End of file 
