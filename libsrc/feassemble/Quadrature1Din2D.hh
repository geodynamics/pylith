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

#include "QuadratureEngine.hh"

namespace pylith {
  namespace feassemble {
    class Quadrature1Din2D;
  } // feassemble
} // pylith

class pylith::feassemble::Quadrature1Din2D : public QuadratureEngine
{ // Quadrature1Din2D
  friend class TestQuadrature1Din2D; // unit testing

// PUBLIC MEMBERS ///////////////////////////////////////////////////////
public :

  /** Constructor.
   *
   * @param q Quadrature information for reference cell.
   */
  Quadrature1Din2D(const QuadratureRefCell& q);

  /// Destructor
  ~Quadrature1Din2D(void);

  /// Create a copy of this object.
  QuadratureEngine* clone(void) const;

  /** Compute geometric quantities for a cell at quadrature points.
   *
   * @param vertCoords Coordinates of vertices of finite-element cell.
   * @param coordDim Spatial dimension of coordinate system.
   * @param cell Finite-element cell
   */
  void computeGeometry(const double* vertCoords,
		       const int coordDim,
		       const int cell);

// PROTECTED METHODS ////////////////////////////////////////////////////
protected :

  /** Copy constructor.
   *
   * @param q Quadrature to copy
   */
  Quadrature1Din2D(const Quadrature1Din2D& q);

// PRIVATE METHODS //////////////////////////////////////////////////////
private :

  /// Not implemented
  const Quadrature1Din2D& operator=(const Quadrature1Din2D&);

}; // Quadrature1Din2D

#include "Quadrature1Din2D.icc" // inline methods

#endif // pylith_feassemble_quadrature1din2d_hh


// End of file 
