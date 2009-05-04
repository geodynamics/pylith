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

#include "QuadratureEngine.hh"

namespace pylith {
  namespace feassemble {
    class Quadrature2D;
    class TestQuadrature2D;
  } // feassemble
} // pylith

class pylith::feassemble::Quadrature2D : public QuadratureEngine
{ // Quadrature2D
  friend class TestQuadrature2D; // unit testing

// PUBLIC MEMBERS ///////////////////////////////////////////////////////
public :

  /** Constructor.
   *
   * @param q Quadrature information for reference cell.
   */
  Quadrature2D(const QuadratureRefCell& q);

  /// Destructor
  ~Quadrature2D(void);

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
  Quadrature2D(const Quadrature2D& q);

// PRIVATE METHODS //////////////////////////////////////////////////////
private :

  const Quadrature2D& operator=(const Quadrature2D&); ///< Not implemented

}; // Quadrature2D

#include "Quadrature2D.icc" // inline methods

#endif // pylith_feassemble_quadrature2d_hh


// End of file 
