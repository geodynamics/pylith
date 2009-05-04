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

#include "QuadratureEngine.hh"

namespace pylith {
  namespace feassemble {
    class Quadrature1Din3D;
  } // feassemble
} // pylith

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
  Quadrature1Din3D(const Quadrature1Din3D& q);

// PRIVATE METHODS //////////////////////////////////////////////////////
private :

  /// Not implemented
  const Quadrature1Din3D& operator=(const Quadrature1Din3D&);

}; // Quadrature1Din3D

#include "Quadrature1Din3D.icc" // inline methods

#endif // pylith_feassemble_quadrature1din3d_hh


// End of file 
