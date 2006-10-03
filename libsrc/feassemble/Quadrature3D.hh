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
 * @file pylith/feassemble/Quadrature3D.hh
 *
 * @brief Quadrature for 3-D finite-elements.
 */

#if !defined(pylith_feassemble_quadrature3d_hh)
#define pylith_feassemble_quadrature3d_hh

#include "Quadrature.hh"

namespace pylith {
  namespace feassemble {
    class Quadrature3D;
    class TestQuadrature3D;
  } // feassemble
} // pylith

class pylith::feassemble::Quadrature3D : public Quadrature
{ // Quadrature1D
  friend class TestQuadrature3D; // unit testing

// PUBLIC MEMBERS ///////////////////////////////////////////////////////
public :

  /// Constructor
  Quadrature3D(void);

  /// Destructor
  virtual ~Quadrature3D(void);

  /// Create a copy of this object.
  virtual
  Quadrature* clone(void) const;

// PROTECTED METHODS ////////////////////////////////////////////////////
protected :

  /** Copy constructor.
   *
   * @param q Quadrature to copy
   */
  Quadrature3D(const Quadrature3D& q);

  /** Compute geometric quantities for a cell.
   *
   * @param coordinates Section containing vertex coordinates
   * @param cell Finite-element cell
   */
  void _computeGeometry(const ALE::Obj<ALE::Mesh::real_section_type>& coordinates,
			const ALE::Mesh::point_type& cell);

// PRIVATE METHODS //////////////////////////////////////////////////////
private :

  /// Not implemented
  const Quadrature3D& operator=(const Quadrature3D&);

}; // Quadrature3D

#include "Quadrature3D.icc" // inline methods

#endif // pylith_feassemble_quadrature3d_hh

// End of file 
