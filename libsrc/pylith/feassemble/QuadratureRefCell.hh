// -*- C++ -*-
//
// ======================================================================
//
// Brad T. Aagaard, U.S. Geological Survey
// Charles A. Williams, GNS Science
// Matthew G. Knepley, University of Chicago
//
// This code was developed as part of the Computational Infrastructure
// for Geodynamics (http://geodynamics.org).
//
// Copyright (c) 2010-2015 University of California, Davis
//
// See COPYING for license information.
//
// ======================================================================
//

/**
 * @file libsrc/feassemble/QuadratureRefCell.hh
 *
 * @brief Object with basic quadrature information for the reference
 * cell.
 */

#if !defined(pylith_feassemble_quadraturerefcell_hh)
#define pylith_feassemble_quadraturerefcell_hh

// Include directives ---------------------------------------------------
#include "feassemblefwd.hh" // forward declarations

#include "CellGeometry.hh" // USES CellGeometry
#include "pylith/utils/array.hh" // HASA scalar_array

// Quadrature -----------------------------------------------------------
/** @brief Object with basic quadrature information for the reference
 * cell.
 *
 * This object contains the basis functions and their derivatives
 * evaluated at the quadrature points of the reference element, and
 * the coordinates and weights of the quadrature points.
 *
 * The Quadrature object manages the general functionality of the
 * numerical quadrature.
 */
class pylith::feassemble::QuadratureRefCell
{ // Quadrature
  friend class TestQuadratureRefCell; // unit testing

// PUBLIC METHODS ///////////////////////////////////////////////////////
public :

  /// Constructor
  QuadratureRefCell(void);

  /// Destructor
  virtual
  ~QuadratureRefCell(void);

  /// Deallocate PETSc and local data structures.
  virtual
  void deallocate(void);
  
  /** Set basis functions and their derivatives, and coordinates and
   *  weights of the quadrature points.
   *
   * @param basis Array of basis functions evaluated at quadrature pts
   *   N0Qp0, N1Qp0, ...
   *   N0Qp1, N1Qp1, ...
   *   ...
   *   size = numQuadPts * numBasis
   *   index = iQuadPt*numBasis + iBasis
   * @param numQuadPts1 Dimension of basis.
   * @param numBasis1 Dimension of basis.
   *
   * @param basisDerivRef Array of basis function derivaties evaluated at
   * quadrature pts, where derivatives are with respect to cell's
   * local coordinates.
   *   N0pQp0, N0qQp0, N0rQp0, N1pQp0, N1qQp0, N1rQp0, ... 
   *   N0pQp1, N0qQp1, N0rQp1, N1pQp1, N1qQp1, N1rQp1, ...
   *   ...
   *   size = numQuadPts * numBasis * cellDim
   *   index = iQuadPt*numBasis*cellDim + iBasis*cellDim + iDim
   * @param numQuadPts2 Dimension of basisDerivRef.
   * @param numBasis2 Dimension of basisDerivRef.
   * @param cellDim2 Dimension of basisDerivRef.
   *
   * @param quadPtsRef Array of coordinates of quadrature points in 
   *   reference cell
   *   Qp0p, Qp0q, Qp0r
   *   Qp1p, Qp1q, Qp1r
   *   size = numQuadPts * cellDim
   *   index = iQuadPt*numDims + iDim
   * @param numQuadPts3 Dimension of quadPtsRef.
   * @param cellDim3 Dimension of quadPtsRef.
   *
   * @param quadWts Array of weights of quadrature points
   *   WtQp0, WtQp1, ...
   *   index = iQuadPt
   * @param numQuadPts4 Dimension of quadWts.
   *
   * @param spaceDim Number of dimensions in coordinates of cell vertices
   */
  void initialize(const PylithScalar* basis,
		  const int numQuadPts1,
		  const int numBasis1,
		  const PylithScalar* basisDerivRef,
		  const int numQuadPts2,
		  const int numBasis2,
		  const int cellDim2,
		  const PylithScalar* quadPtsRef,
		  const int numQuadPts3,
		  const int cellDim3,
		  const PylithScalar* quadWts,
		  const int numQuadPts4,
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
  void minJacobian(const PylithScalar min);

  /** Get minimum allowable determinant of Jacobian.
   *
   * @returns Minimum allowable value for Jacobian
   */
  PylithScalar minJacobian(void) const;

  /** Get coordinates of quadrature points in reference cell.
   *
   * @returns Array of coordinates of quadrature points in reference cell.
   */
  const scalar_array& quadPtsRef(void) const;

  /** Get weights of quadrature points.
   *
   * @returns Weights of quadrature points
   */
  const scalar_array& quadWts(void) const;

  /** Get basis fns evaluated at quadrature points.
   *
   * @returns Array of basis fns evaluated at quadrature points
   */
  const scalar_array& basis(void) const;

  /** Get derivates of basis fns evaluated at quadrature points.
   *
   * @returns Array of derivates of basis fns evaluated at quadrature points
   */
  const scalar_array& basisDerivRef(void) const;

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

  /** Compute minimum width across cell.
   *
   * @param coordinatesCell Coordinates of vertices in cell.
   * @param numVertices Number of vertices in cell.
   * @param spaceDim Coordinate dimension.
   *
   * @returns Minimum width across cell.
   */
  PylithScalar minCellWidth(const PylithScalar* coordinatesCell,
			    const int numVertices,
			    const int spaceDim) const;

// PROTECTED METHODS ////////////////////////////////////////////////////
protected :

  /** Copy constructor.
   *
   * @param q Quadrature to copy
   */
  QuadratureRefCell(const QuadratureRefCell& q);

// PROTECTED MEMBERS ////////////////////////////////////////////////////
protected :

  PylithScalar _minJacobian; ///< Minium allowable Jacobian determinant
  
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
  scalar_array _quadPtsRef;

  /** Array of weights of quadrature points.
   *
   * WtQp0, WtQp1, ...
   * size = numQuadPts
   * index = iQuadPt
   */
  scalar_array _quadWts;

  /** Array of basis functions evaluated at the quadrature points.
   *
   * N0Qp0, N1Qp0, ...
   * N0Qp1, N1Qp1, ...
   *
   * size = numQuadPts * numBasis
   * index = iQuadPt*numBasis + iBasis
   */
  scalar_array _basis;

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
  scalar_array _basisDerivRef;

  int _cellDim; ///< Number of dimensions in reference cell
  int _numBasis; ///< Number of basis functions (and vertices) for cell
  int _numQuadPts; ///< Number of quadrature points
  int _spaceDim; ///< Number of dimensions in coordinates of cell vertices

  CellGeometry* _geometry; ///< Geometry of reference cell

// NOT IMPLEMENTED //////////////////////////////////////////////////////
private :

  const QuadratureRefCell& operator=(const QuadratureRefCell&); ///< Not implemented

}; // QuadratureRefCell

#include "QuadratureRefCell.icc" // inline methods

#endif // pylith_feassemble_quadraturerefcell_hh


// End of file 
