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
 * @file pylith/feassemble/Quadrature0D.hh
 *
 * @brief Quadrature for 0-D finite-elements.
 *
 * Need Quadrature in 0-D for integration of boundary condition for
 * 1-D meshes.
 */

#if !defined(pylith_feassemble_quadrature0d_hh)
#define pylith_feassemble_quadrature0d_hh

#include "QuadratureEngine.hh"

namespace pylith {
  namespace feassemble {
    class TestQuadrature0D;
  } // feassemble
} // pylith

class pylith::feassemble::Quadrature0D : public QuadratureEngine
{ // Quadrature0D
  friend class TestQuadrature0D; // unit testing

// PUBLIC MEMBERS ///////////////////////////////////////////////////////
public :

  /** Constructor.
   *
   * @param q Quadrature information for reference cell.
   */
  Quadrature0D(const QuadratureRefCell& q);

  /// Destructor
  ~Quadrature0D(void);

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
  Quadrature0D(const Quadrature0D& q);

// PRIVATE METHODS //////////////////////////////////////////////////////
private :

  const Quadrature0D& operator=(const Quadrature0D&); ///< Not implemented

}; // Quadrature0D

#include "Quadrature0D.icc" // inline methods

#endif // pylith_feassemble_quadrature0d_hh


// End of file 
