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
 * @file pylith/feassemble/Quadrature1Din3D.hh
 *
 * @brief Quadrature for 1-D finite-elements in 3-D space.
 */

#if !defined(pylith_feassemble_quadrature1din3d_hh)
#define pylith_feassemble_quadrature1din3d_hh

#include "Quadrature.hh"

namespace pylith {
  namespace feassemble {
    class Quadrature1Din3D;
    class TestQuadrature1Din3D;
  } // feassemble
} // pylith

class pylith::feassemble::Quadrature1Din3D : public Quadrature
{ // Quadrature1D
  friend class TestQuadrature1Din3D; // unit testing

// PUBLIC MEMBERS ///////////////////////////////////////////////////////
public :

  /// Constructor
  Quadrature1Din3D(void);

  /// Destructor
  virtual ~Quadrature1Din3D(void);

  /// Create a copy of this object.
  virtual
  Quadrature* clone(void) const;

// PROTECTED METHODS ////////////////////////////////////////////////////
protected :

  /** Copy constructor.
   *
   * @param q Quadrature to copy
   */
  Quadrature1Din3D(const Quadrature1Din3D& q);

  /** Compute geometric quantities for a cell.
   *
   * @param coordinates Section containing vertex coordinates
   * @param cell Finite-element cell
   */
  void _computeGeometry(const ALE::Obj<ALE::Mesh::section_type>& coordinates,
			const ALE::Mesh::point_type& cell);

// PRIVATE METHODS //////////////////////////////////////////////////////
private :

  /// Not implemented
  const Quadrature1Din3D& operator=(const Quadrature1Din3D&);

}; // Quadrature1Din3D

#include "Quadrature1Din3D.icc" // inline methods

#endif // pylith_feassemble_quadrature1din3d_hh

// End of file 
