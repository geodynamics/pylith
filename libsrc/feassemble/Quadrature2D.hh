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
{ // Quadrature2D
  friend class TestQuadrature2D; // unit testing

// PUBLIC MEMBERS ///////////////////////////////////////////////////////
public :

  /// Constructor
  Quadrature2D(void);

  /// Destructor
  ~Quadrature2D(void);

  /// Create a copy of this object.
  Quadrature* clone(void) const;

  /** Compute geometric quantities for a cell.
   *
   * @param coordinates Section containing vertex coordinates
   * @param cell Finite-element cell
   */
  void computeGeometry(const ALE::Obj<real_section_type>& coordinates,
		       const topology_type::point_type& cell);

// PROTECTED METHODS ////////////////////////////////////////////////////
protected :

  /** Copy constructor.
   *
   * @param q Quadrature to copy
   */
  Quadrature2D(const Quadrature2D& q);

// PRIVATE METHODS //////////////////////////////////////////////////////
private :

  const Quadrature2D& operator=(const Quadrature2D&); ///< Not implemented

}; // Quadrature2D

#include "Quadrature2D.icc" // inline methods

#endif // pylith_feassemble_quadrature2d_hh

// End of file 
