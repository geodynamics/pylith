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

#include "pylith/utils/array.hh" // HASA double_array
#include "pylith/utils/sievetypes.hh" // USES PETSc Mesh

namespace pylith {
  namespace feassemble {
    class Quadrature;

    class CellGeometry; // HOLDSA CellGeometry
    class TestQuadrature;
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

  /** Set vertices (dual basis), basis functions and their
   *  derivatives, and coordinates and weights of the quadrature
   *  points.
   *
   * @param vertices Array of coordinates of vertices (dual basis).
   *   v0x, v0y, v0z
   *   v1x, v1y, v1z
   *   ...
   *   size = numBasis * cellDim
   *   index = iBasis*cellDim + iDim
   *
   * @param basis Array of basis functions evaluated at quadrature pts
   *   N0Qp0, N1Qp0, ...
   *   N0Qp1, N1Qp1, ...
   *   ...
   *   size = numQuadPts * numBasis
   *   index = iQuadPt*numBasis + iBasis
   *
   * @param basisDerivRef Array of basis function derivaties evaluated at
   * quadrature pts, where derivatives are with respect to cell's
   * local coordinates.
   *   N0pQp0, N0qQp0, N0rQp0, N1pQp0, N1qQp0, N1rQp0, ... 
   *   N0pQp1, N0qQp1, N0rQp1, N1pQp1, N1qQp1, N1rQp1, ...
   *   ...
   *   size = numQuadPts * numBasis * cellDim
   *   index = iQuadPt*numBasis*cellDim + iBasis*cellDim + iDim
   *
   * @param quadPts Array of coordinates of quadrature points in 
   *   reference cell
   *   Qp0p, Qp0q, Qp0r
   *   Qp1p, Qp1q, Qp1r
   *   size = numQuadPts * numDims
   *   index = iQuadPt*numDims + iDim
   *
   * @param quadWts Array of weights of quadrature points
   *   WtQp0, WtQp1, ...
   *   index = iQuadPt
   *
   * @param cellDim Number of dimensions in reference cell
   * @param numBasis Number of basis functions (and vertices) for a cell
   * @param numQuadPts Number of quadrature points
   * @param spaceDim Number of dimensions in coordinates of cell vertices
   */
  void initialize(const double* vertices,
		  const double* basis,
		  const double* basisDerivRef,
		  const double* quadPtsRef,
		  const double* quadWts,
		  const int cellDim,
		  const int numBasis,
		  const int numQuadPts,
		  const int spaceDim);

  /** Set geometry associated with reference cell.
   *
   * @param geometry Geometry of reference cell.
   */
  void refGeometry(CellGeometry* const geometry);

  /** Get geometry associated with reference cell.
   *
   * @returns Geometry of reference cell.
   */
  const CellGeometry& refGeometry(void) const;

  /** Set minimum allowable determinant of Jacobian.
   *
   * @param tolerance Minimum allowable value for Jacobian
   */
  void minJacobian(const double min);

  /** Get minimum allowable determinant of Jacobian.
   *
   * @returns Minimum allowable value for Jacobian
   */
  double minJacobian(void);

  /** Get coordinates of vertices in reference cell (dual basis).
   *
   * @returns Array of coordinates of vertices in reference cell
   */
  const double_array& vertices(void) const;

  /** Get coordinates of quadrature points in cell (NOT reference cell).
   *
   * @returns Array of coordinates of quadrature points in cell
   */
  const double_array& quadPts(void) const;

  /** Get weights of quadrature points.
   *
   * @returns Weights of quadrature points
   */
  const double_array& quadWts(void) const;

  /** Get basis fns evaluated at quadrature points.
   *
   * @returns Array of basis fns evaluated at quadrature points
   */
  const double_array& basis(void) const;

  /** Get derivatives of basis fns evaluated at quadrature points.
   *
   * @returns Array of derivatives of basis fns evaluated at
   * quadrature points
   */
  const double_array& basisDeriv(void) const;

  /** Get Jacobians evaluated at quadrature points.
   *
   * @returns Array of Jacobian inverses evaluated at quadrature points.
   */
  const double_array& jacobian(void) const;

  /** Get determinants of Jacobian evaluated at quadrature points.
   *
   * @returns Array of determinants of Jacobian evaluated at quadrature pts
   */
  const double_array& jacobianDet(void) const;

  /** Get Jacobian inverses evaluated at quadrature points.
   *
   * @returns Array of Jacobian inverses evaluated at quadrature points.
   */
  const double_array& jacobianInv(void) const;

  /** Get number of dimensions in reference cell.
   *
   * @returns Number of dimensions in reference cell
   */
  int cellDim(void) const;

  /** Get number of basis functions for cell.
   *
   * @returns Number of basis functions for cell
   */
  int numBasis(void) const;

  /** Get number of quadrature points.
   *
   * @returns Number of quadrature points
   */
  int numQuadPts(void) const;

  /** Get number of dimensions in coordinates of cell vertices.
   *
   * @returns Number of dimensions in coordinates of cell vertices
   */
  int spaceDim(void) const;

  /** Compute geometric quantities for a cell at quadrature points.
   *
   * @param mesh Finite-element mesh
   * @param coordinates Section containing vertex coordinates
   * @param cell Finite-element cell
   */
  virtual 
  void computeGeometry(const ALE::Obj<Mesh>& mesh,
		       const ALE::Obj<real_section_type>& coordinates,
		       const Mesh::point_type& cell) = 0;

// PROTECTED METHODS ////////////////////////////////////////////////////
protected :

  /** Copy constructor.
   *
   * @param q Quadrature to copy
   */
  Quadrature(const Quadrature& q);

  /// Check determinant of Jacobian against minimum allowable value
  void _checkJacobianDet(const double det) const;

  /// Set entries in geometry arrays to zero.
  void _resetGeometry(void);

// PRIVATE METHODS //////////////////////////////////////////////////////
private :

  const Quadrature& operator=(const Quadrature&); ///< Not implemented

// PROTECTED MEMBERS ////////////////////////////////////////////////////
protected :

  double _minJacobian; ///< Minium allowable Jacobian determinant
  
  /** Array of coordinates of vertices in reference cell (dual basis).
   *
   * Reference coordinates: (p,q,r)
   *
   * v0p, v0q, v0r
   * v1p, v1q, v1r
   *
   * size = numBasis * cellDim
   * index = iBasis*cellDim + iDim
   */
  double_array _vertices;

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
  double_array _quadPtsRef;

  /** Array of coordinates of quadrature points in cell (global coordinates).
   *
   * Qp0x, Qp0y, Qp0z
   * Qp1x, Qp1y, Qp1z
   *
   * size = numQuadPts * spaceDim
   * index = iQuadPts*spaceDim + iDim
   */
  double_array _quadPts;

  /** Array of weights of quadrature points.
   *
   * WtQp0, WtQp1, ...
   * size = numQuadPts
   * index = iQuadPt
   */
  double_array _quadWts;

  /** Array of basis functions evaluated at the quadrature points.
   *
   * N0Qp0, N1Qp0, ...
   * N0Qp1, N1Qp1, ...
   *
   * size = numQuadPts * numBasis
   * index = iQuadPt*numBasis + iBasis
   */
  double_array _basis;

  /** Array of basis function derivatives evaluated at the quadrature
   * points, where derivatives are with respect to cell's local
   * coordinates.
   *
   * N0pQp0, N0qQp0, N0rQp0, N1pQp0, N1qQp0, N1rQp0, ... 
   * N0pQp1, N0qQp1, N0rQp1, N1pQp1, N1qQp1, N1rQp1, ...
   *
   * size = numQuadPts * numBasis * spaceDim
   * index = iQuadPt*numBasis*cellDim + iBasis*spaceDim + iDim
   */
  double_array _basisDerivRef;

  /** Array of basis function derivatives evaluated at the quadrature
   * points, where derivatives are with respect to global coordinates.
   *
   * N0xQp0, N0yQp0, N0zQp0, N1xQp0, N1yQp0, N1zQp0, ... 
   * N0xQp1, N0yQp1, N0zQp1, N1xQp1, N1yQp1, N1zQp1, ...
   *
   * size = numQuadPts * numBasis * cellDim
   * index = iQuadPt*numBasis*cellDim + iBasis*cellDim + iDim
   */
  //double_array _basisDeriv;

  /** Array of Jacobians evaluated at quadrature points.
   *
   * Qp0J00, Qp0J01, Qp0J02, ...
   * Qp1J00, Qp1J01, Qp1J02, ...
   * ...
   *
   * size = numQuadPts*cellDim*spaceDim
   * index = iQuadPt*cellDim*spaceDim + iRow*spaceDim + iCol
   */
  double_array _jacobian;

  /** Array of determinant of Jacobian evaluated at quadrature points.
   *
   * JdetQp0, JdetQp1, ...
   *
   * size = numQuadPts
   * index = iQuadPt
   */
  double_array _jacobianDet;

  /** Array of Jacobian inverses evaluated at quadrature points.
   *
   * Qp0Jinv00, Qp0Jinv01, Qp0Jinv02, ...
   * Qp1Jinv00, Qp1Jinv01, Qp1Jinv02, ...
   * ...
   *
   * size = numQuadPts*spaceDim*cellDim
   * index = iQuadPt*spaceDim*cellDim + iRow*cellDim + iCol
   */
  double_array _jacobianInv;

  int _cellDim; ///< Number of dimensions in reference cell
  int _numBasis; ///< Number of basis functions (and vertices) for cell
  int _numQuadPts; ///< Number of quadrature points
  int _spaceDim; ///< Number of dimensions in coordinates of cell vertices

  CellGeometry* _geometry; ///< Geometry of reference cell

}; // Quadrature

#include "Quadrature.icc" // inline methods

#endif // pylith_feassemble_quadrature_hh

// End of file 
