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
 * @file pylith/feassemble/Quadrature0D.hh
 *
 * @brief Quadrature for 0-D finite-elements.
 *
 * Need Quadrature in 0-D for integration of boundary condition for
 * 1-D meshes.
 */

#if !defined(pylith_feassemble_quadrature0d_hh)
#define pylith_feassemble_quadrature0d_hh

#include "Quadrature.hh"

namespace pylith {
  namespace feassemble {
    class Quadrature0D;
    class TestQuadrature0D;
  } // feassemble
} // pylith

class pylith::feassemble::Quadrature0D : public Quadrature
{ // Quadrature0D
  friend class TestQuadrature0D; // unit testing

// PUBLIC MEMBERS ///////////////////////////////////////////////////////
public :

  /// Constructor
  Quadrature0D(void);

  /// Destructor
  ~Quadrature0D(void);

  /// Create a copy of this object.
  Quadrature* clone(void) const;

  /** Compute geometric quantities for a cell at quadrature points.
   *
   * @param mesh Finite-element mesh
   * @param coordinates Section containing vertex coordinates
   * @param cell Finite-element cell
   */
  void computeGeometry(const ALE::Obj<Mesh>& mesh,
		       const ALE::Obj<real_section_type>& coordinates,
		       const Mesh::point_type& cell);

// PROTECTED METHODS ////////////////////////////////////////////////////
protected :

  /** Copy constructor.
   *
   * @param q Quadrature to copy
   */
  Quadrature0D(const Quadrature0D& q);

// PRIVATE METHODS //////////////////////////////////////////////////////
private :

  const Quadrature0D& operator=(const Quadrature0D&); ///< Not implemented

}; // Quadrature0D

#include "Quadrature0D.icc" // inline methods

#endif // pylith_feassemble_quadrature0d_hh

// End of file 
