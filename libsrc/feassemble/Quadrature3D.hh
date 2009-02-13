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

#include "QuadratureEngine.hh" // ISA QuadratureEngine

class pylith::feassemble::Quadrature3D : public QuadratureEngine
{ // Quadrature3D
  friend class TestQuadrature3D; // unit testing

// PUBLIC MEMBERS ///////////////////////////////////////////////////////
public :

  /** Constructor.
   *
   * @param q Quadrature information for reference cell.
   */
  Quadrature3D(const QuadratureRefCell& q);

  /// Destructor
  ~Quadrature3D(void);

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
  Quadrature3D(const Quadrature3D& q);

// PRIVATE METHODS //////////////////////////////////////////////////////
private :

  /// Not implemented
  const Quadrature3D& operator=(const Quadrature3D&);

}; // Quadrature3D

#include "Quadrature3D.icc" // inline methods

#endif // pylith_feassemble_quadrature3d_hh


// End of file 
