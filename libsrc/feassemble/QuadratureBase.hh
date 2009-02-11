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
 * @file pylith/feassemble/QuadratureBase.hh
 *
 * @brief Object with basic quadrature information for the reference cell.
 *
 * This object contains the basis functions and their derivatives
 * evaluated at the quadrature points of the reference element, and
 * the coordinates and weights of the quadrature points. 
 *
 * The Quadrature object manages the general functionality of the
 * numerical quadrature.
 */

#if !defined(pylith_feassemble_quadraturebase_hh)
#define pylith_feassemble_quadraturebase_hh

// Include directives ---------------------------------------------------
#include "feassemblefwd.hh" // forward declarations

#include "pylith/utils/array.hh" // HASA double_array

// Quadrature -----------------------------------------------------------
class pylith::feassemble::QuadratureBase
{ // Quadrature
  friend class TestQuadratureBase; // unit testing

// PUBLIC METHODS ///////////////////////////////////////////////////////
public :

  /// Constructor
  QuadratureBase(void);

  /// Destructor
  ~QuadratureBase(void);

  /** Set basis functions and their derivatives, and coordinates and
   *  weights of the quadrature points.
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
   * @param numBasis Number of basis functions for a cell
   * @param numQuadPts Number of quadrature points
   * @param spaceDim Number of dimensions in coordinates of cell vertices
   */
  void initialize(const double* basis,
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
  double minJacobian(void) const;

  /** Get coordinates of quadrature points in reference cell.
   *
   * @returns Array of coordinates of quadrature points in reference cell.
   */
  const double_array& quadPtsRef(void) const;

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

// PROTECTED METHODS ////////////////////////////////////////////////////
protected :

  /** Copy constructor.
   *
   * @param q Quadrature to copy
   */
  QuadratureBase(const QuadratureBase& q);

  /* Check determinant of Jacobian against minimum allowable value.
   *
   * @param det Value of determinant of Jacobian
   * @param cell Label of finite-element cell
   */
  void _checkJacobianDet(const double det,
			 const int cell) const;

// PROTECTED MEMBERS ////////////////////////////////////////////////////
protected :

  double _minJacobian; ///< Minium allowable Jacobian determinant
  
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
   * size = numQuadPts * numBasis * cellDim
   * index = iQuadPt*numBasis*cellDim + iBasis*cellDim + iDim
   */
  double_array _basisDerivRef;

  int _cellDim; ///< Number of dimensions in reference cell
  int _numBasis; ///< Number of basis functions (and vertices) for cell
  int _numQuadPts; ///< Number of quadrature points
  int _spaceDim; ///< Number of dimensions in coordinates of cell vertices

  CellGeometry* _geometry; ///< Geometry of reference cell

// NOT IMPLEMENTED //////////////////////////////////////////////////////
private :

  const QuadratureBase& operator=(const QuadratureBase&); ///< Not implemented

}; // QuadratureBase

#include "QuadratureBase.icc" // inline methods

#endif // pylith_feassemble_quadraturebase_hh


// End of file 
