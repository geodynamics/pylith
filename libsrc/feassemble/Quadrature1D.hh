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
 * @file pylith/feassemble/Quadrature1D.hh
 *
 * @brief Quadrature for 1-D finite-elements.
 */

#if !defined(pylith_feassemble_quadrature1d_hh)
#define pylith_feassemble_quadrature1d_hh

#include "Quadrature.hh"

namespace pylith {
  namespace feassemble {
    class Quadrature1D;
    class TestQuadrature1D;
  } // feassemble
} // pylith

class pylith::feassemble::Quadrature1D : public Quadrature
{ // Quadrature1D
  friend class TestQuadrature1D; // unit testing

// PUBLIC MEMBERS ///////////////////////////////////////////////////////
public :

  /// Constructor
  Quadrature1D(void);

  /// Destructor
  virtual ~Quadrature1D(void);

  /// Create a copy of this object.
  virtual
  Quadrature* clone(void) const;

// PROTECTED METHODS ////////////////////////////////////////////////////
protected :

  /** Copy constructor.
   *
   * @param q Quadrature to copy
   */
  Quadrature1D(const Quadrature1D& q);

  /** Compute geometric quantities for a cell.
   *
   * @param coordinates Section containing vertex coordinates
   * @param cell Finite-element cell
   */
  void _computeGeometry(const ALE::Obj<ALE::Mesh::real_section_type>& coordinates,
			const ALE::Mesh::point_type& cell);

// PRIVATE METHODS //////////////////////////////////////////////////////
private :

  const Quadrature1D& operator=(const Quadrature1D&); ///< Not implemented

}; // Quadrature1D

#include "Quadrature1D.icc" // inline methods

#endif // pylith_feassemble_quadrature1d_hh

// End of file 
