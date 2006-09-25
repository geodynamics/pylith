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
 * @file pylith/feassemble/Quadrature1Din2D.hh
 *
 * @brief Quadrature for 1-D finite-elements in 2-D space.
 */

#if !defined(pylith_feassemble_quadrature1din2d_hh)
#define pylith_feassemble_quadrature1din2d_hh

#include "Quadrature.hh"

namespace pylith {
  namespace feassemble {
    class Quadrature1Din2D;
    class TestQuadrature1Din2D;
  } // feassemble
} // pylith

class pylith::feassemble::Quadrature1Din2D : public Quadrature
{ // Quadrature1D
  friend class TestQuadrature1Din2D; // unit testing

// PUBLIC MEMBERS ///////////////////////////////////////////////////////
public :

  /// Constructor
  Quadrature1Din2D(void);

  /// Destructor
  virtual ~Quadrature1Din2D(void);

  /// Create a copy of this object.
  virtual
  Quadrature* clone(void) const;

// PROTECTED METHODS ////////////////////////////////////////////////////
protected :

  /** Copy constructor.
   *
   * @param q Quadrature to copy
   */
  Quadrature1Din2D(const Quadrature1Din2D& q);

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
  const Quadrature1Din2D& operator=(const Quadrature1Din2D&);

}; // Quadrature1Din2D

#include "Quadrature1Din2D.icc" // inline methods

#endif // pylith_feassemble_quadrature1din2d_hh

// End of file 
