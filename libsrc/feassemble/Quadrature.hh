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
 * @file pylith/feassemble/Quadrature.hh
 *
 * @brief Abstract base class for integrating over finite-elements
 * using quadrature.
 *
 * This object contains the basis functions and their derivatives
 * evaluated at the quadrature points of the reference element, and
 * the coordinates and weights of the quadrature points. Given a cell
 * this object will compute the cell's Jacobian, the determinant of
 * the Jacobian, the inverse of the Jacobian, and the coordinates in
 * the domain of the cell's quadrature points. The Jacobian and its
 * inverse are computed at the quadrature points.
 *
 * The memory for the Jacobian and its associated information are
 * managed locally.
 */

#if !defined(pylith_feassemble_quadrature_hh)
#define pylith_feassemble_quadrature_hh

#include <Mesh.hh>

namespace pylith {
  namespace feassemble {
    class Quadrature;
  } // feassemble
} // pylith

class pylith::feassemble::Quadrature
{ // Quadrature
  friend class TestQuadrature; // unit testing

// PUBLIC METHODS ///////////////////////////////////////////////////////
public :

  /// Constructor
  Quadrature(void);

  /// Destructor
  virtual
  ~Quadrature(void);

  /// Create a copy of this object.
  virtual
  Quadrature* clone(void) const = 0;

  /** Set basis functions and their derivatives, and coordinates and
   *  weights of the quadrature points.
   *
   * @param basis Array of basis functions evaluated at quadrature pts
   *   N0Qp0, N1Qp0, ...
   *   N0Qp1, N1Qp1, ...
   *   ...
   *   size = numQuadPts * numCorners
   *   index = iQuadPt*numCorners + iBasis
   *
   * @param basisDeriv Array of basis function derivaties evaluated 
   *   at quadrature pts
   *   N0xQp0, N0yQp0, N0zQp0, N1xQp0, N1yQp0, N1zQp0, ... 
   *   N0xQp1, N0yQp1, N0zQp1, N1xQp1, N1yQp1, N1zQp1, ...
   *   ...
   *   size = numCorners * numQuadPts * spaceDim
   *   index = iQuadPt*numCorners*spaceDim + iBasis*spaceDim + iDim
   *
   * @param quadPts Array of coordinates of quadrature points in 
   *   reference cell
   *   Qp0x, Qp0y, Qp0z
   *   Qp1x, Qp1y, Qp1z
   *   size = numQuadPts * numDims
   *   index = iQuadPt*numDims + iDim
   *
   * @param quadWts Array of weights of quadrature points
   *   WtQp0, WtQp1, ...
   *   index = iQuadPt
   *
   * @param cellDim Number of dimensions in reference cell
   * @param numCorners Number of vertices in a cell
   * @param numQuadPts Number of quadrature points
   * @param spaceDim Number of dimensions in coordinates of cell vertices
   */
  void initialize(const double* basis,
		  const double* basisDeriv,
		  const double* quadPtsRef,
		  const double* quadWts,
		  const int cellDim,
		  const int numCorners,
		  const int numQuadPts,
		  const int spaceDim);

  /** Set tolerance for minimum allowable Jacobian.
   *
   * @param tolerance Minimum allowable value for Jacobian
   */
  void jacobianTolerance(const double tolerance);

  /** Get tolerance for minimum allowable Jacobian.
   *
   * @returns Minimum allowable value for Jacobian
   */
  double jacobianTolerance(void);

// PROTECTED METHODS ////////////////////////////////////////////////////
protected :

  /** Copy constructor.
   *
   * @param q Quadrature to copy
   */
  Quadrature(const Quadrature& q);

  /** Compute geometric quantities for a cell.
   *
   * @param coordinates Section containing vertex coordinates
   * @param cell Finite-element cell
   */
  virtual 
  void _computeGeometry(const ALE::Obj<ALE::Mesh::section_type>& coordinates,
			const ALE::Mesh::point_type& cell) = 0;

// PRIVATE METHODS //////////////////////////////////////////////////////
private :

  const Quadrature& operator=(const Quadrature&); ///< Not implemented

// PROTECTED MEMBERS ////////////////////////////////////////////////////
protected :

  double _jacobianTol; ///< Tolerance for minium allowable Jacobian determinant
  
  /** Array of basis functions evaluated at the quadrature points.
   *
   * N0Qp0, N1Qp0, ...
   * N0Qp1, N1Qp1, ...
   *
   * size = numQuadPts * numCorners
   * index = iQuadPt*numCorners + iBasis
   */
  double* _basis; ///< Array of basis fns evaluated at quad pts

  /** Array of basis functions evaluated at the quadrature points.
   *
   * N0xQp0, N0yQp0, N0zQp0, N1xQp0, N1yQp0, N1zQp0, ... 
   * N0xQp1, N0yQp1, N0zQp1, N1xQp1, N1yQp1, N1zQp1, ...
   *
   * size = numQuadPts * numCorners * spaceDim
   * index = iQuadPt*numCorners*spaceDim + iBasis*spaceDim + iDim
   */
  double* _basisDeriv;

  /** Array of coordinates of quadrature points in reference cell.
   *
   * Reference coordinates: (p,q,r)
   *
   * Qp0p, Qp0q, Qp0r
   * Qp1p, Qp1q, Qp1r
   *
   * size = numQuadPts * cellDim
   * index = iQuadPts*cellDim + iDim
   */
  double* _quadPtsRef;

  /** Array of coordinates of quadrature points in cell (NOT reference cell).
   *
   * Qp0x, Qp0y, Qp0z
   * Qp1x, Qp1y, Qp1z
   *
   * size = numQuadPts * spaceDim
   * index = iQuadPts*spaceDim + iDim
   */
  double* _quadPts;

  /** Array of weights of quadrature points.
   *
   * WtQp0, WtQp1, ...
   * size = numQuadPts
   * index = iQuadPt
   */
  double* _quadWts;

  /** Array of Jacobian evaluated at quadrature points.
   *
   * Qp0J00, Qp0J01, Qp0J02, ...
   * Qp1J00, Qp1J01, Qp1J02, ...
   * ...
   *
   * size = numQuadPts*cellDim*spaceDim
   * index = iQuadPt*cellDim*spaceDim + iRow*spaceDim + iCol
   */
  double* _jacobian;

  /** Array of Jacobian inverses evaluated at quadrature points.
   *
   * Qp0Jinv00, Qp0Jinv01, Qp0Jinv02, ...
   * Qp1Jinv00, Qp1Jinv01, Qp1Jinv02, ...
   * ...
   *
   * size = numQuadPts*spaceDim*cellDim
   * index = iQuadPt*spaceDim*cellDim + iRow*cellDim + iCol
   */
  double* _jacobianInv;

  /** Array of determinant of Jacobian evaluated at quadrature points.
   *
   * JdetQp0, JdetQp1, ...
   *
   * size = numQuadPts
   * index = iQuadPt
   */
  double* _jacobianDet;

  int _cellDim; ///< Number of dimensions in reference cell
  int _numCorners; ///< Number of vertices in cell
  int _numQuadPts; ///< Number of quadrature points
  int _spaceDim; ///< Number of dimensions in coordinates of cell vertices

}; // Quadrature

#include "Quadrature.icc" // inline methods

#endif // pylith_feassemble_quadrature_hh

// End of file 
