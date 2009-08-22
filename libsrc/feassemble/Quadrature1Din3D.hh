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
 * @file libsrc/feassemble/Quadrature1Din3D.hh
 *
 * @brief Quadrature for 1-D finite-elements in 3-D space.
 */

#if !defined(pylith_feassemble_quadrature1din3d_hh)
#define pylith_feassemble_quadrature1din3d_hh

// Include directives ---------------------------------------------------
#include "QuadratureEngine.hh"

// Quadrature0D ---------------------------------------------------------
/** @brief Quadrature for 1-D finite-elements in 3-D space.
 */
class pylith::feassemble::Quadrature1Din3D : public QuadratureEngine
{ // Quadrature1Din3D
  friend class TestQuadrature1Din3D; // unit testing

// PUBLIC MEMBERS ///////////////////////////////////////////////////////
public :

  /** Constructor.
   *
   * @param q Quadrature information for reference cell.
   */
  Quadrature1Din3D(const QuadratureRefCell& q);

  /// Destructor
  ~Quadrature1Din3D(void);

  /// Create a copy of this object.
  QuadratureEngine* clone(void) const;

  /** Compute geometric quantities for a cell at quadrature points.
   *
   * @param coordinatesCell Coordinates of cell's vertices.
   * @param cell Finite-element cell
   */
  void computeGeometry(const double_array& coordinatesCell,
		       const int cell);

// PROTECTED METHODS ////////////////////////////////////////////////////
protected :

  /** Copy constructor.
   *
   * @param q Quadrature to copy
   */
  Quadrature1Din3D(const Quadrature1Din3D& q);

// PRIVATE METHODS //////////////////////////////////////////////////////
private :

  /// Not implemented
  const Quadrature1Din3D& operator=(const Quadrature1Din3D&);

}; // Quadrature1Din3D

#include "Quadrature1Din3D.icc" // inline methods

#endif // pylith_feassemble_quadrature1din3d_hh


// End of file 
