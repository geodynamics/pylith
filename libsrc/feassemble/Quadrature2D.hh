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
 * @file pylith/feassemble/Quadrature2D.hh
 *
 * @brief Quadrature for 2-D finite-elements in 2-D space.
 */

#if !defined(pylith_feassemble_quadrature2d_hh)
#define pylith_feassemble_quadrature2d_hh

#include "Quadrature.hh"

namespace pylith {
  namespace feassemble {
    class Quadrature2D;
    class TestQuadrature2D;
  } // feassemble
} // pylith

class pylith::feassemble::Quadrature2D : public Quadrature
{ // Quadrature1D
  friend class TestQuadrature2D; // unit testing

// PUBLIC MEMBERS ///////////////////////////////////////////////////////
public :

  /// Constructor
  Quadrature2D(void);

  /// Destructor
  virtual ~Quadrature2D(void);

  /// Create a copy of this object.
  virtual
  Quadrature* clone(void) const;

// PROTECTED METHODS ////////////////////////////////////////////////////
protected :

  /** Copy constructor.
   *
   * @param q Quadrature to copy
   */
  Quadrature2D(const Quadrature2D& q);

  /** Compute geometric quantities for a cell.
   *
   * @param coordinates Section containing vertex coordinates
   * @param cell Finite-element cell
   */
  void _computeGeometry(const ALE::Obj<ALE::Mesh::section_type>& coordinates,
			const ALE::Mesh::point_type& cell);

// PRIVATE METHODS //////////////////////////////////////////////////////
private :

  const Quadrature2D& operator=(const Quadrature2D&); ///< Not implemented

}; // Quadrature2D

#include "Quadrature2D.icc" // inline methods

#endif // pylith_feassemble_quadrature2d_hh

// End of file 
