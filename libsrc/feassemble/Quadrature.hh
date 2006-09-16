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

// @file pylith/feassemble/Quadrature.hh

// @brief Abstract base class for integrating over finite-elements
// using quadratures.
//
// This object holds the basis functions and their derivatives
// evaluated at the quadrature points, and the coordinates and weights
// of the quadrature points.

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
  
// PUBLIC MEMBERS ///////////////////////////////////////////////////////
public :

  /// Constructor
  Quadrature(void);

  /// Destructor
  virtual ~Quadrature(void);

  /** Compute geometric quantities for a cell.
   *
   * @param coordinates Section containing vertex coordinates
   * @param cell Finite-element cell
   * @param pV Array of ???
   * @param pJacobian Jacobian evaluated at quadrature points
   *   size = numDims*numDims*numQuadPts
   *   index = iQuadPt*numDims*numDims + iJacobian
   * @param pJacobianInv Inverse Jacobian evaluated at quadrature points
   *   quadrature pts
   *   size = numDims*numDims*numQuadPts
   *   index = iQuadPt*numDims*numDims + iJacobian
   * @param jacobianDet Determinant of Jacobian
   */
  virtual void compute(const ALE::Obj<ALE::Mesh::section_type>& coordinates,
		       const ALE::Mesh::point_type& cell,
		       double* pV,
		       double* pJacobian,
		       double* pJacobianInv,
		       double& jacobianDet) = 0;

  /** Set basis functions and their derivatives and coordinates and
   *  weights of the quadrature points.
   *
   * @param pBasisFns Array of basis functions evaluated at quadrature pts
   *   index = iVertex*numDims + iDimension
   *
   * @param pBasisFnsDeriv Array of basis function derivaties evaluated 
   *   at quadrature pts
   *   index = iVertex*numDims + iDimension??
   *
   * @param pQuadPts Array of coordinates of quadrature points in 
   *   reference element
   *   index = iQuadPt*numDims + iDimension
   *
   * @param pQuadWts Array of weights of quadrature points
   *   index = iQuadPt
   *
   * @param numDims Number of dimensions
   * @param numCorners Number of vertices in a cell
   * @param numQuadPts Number of quadrature points
   */
  void initialize(const double* pBasisFns,
		  const double* pBasisFnsDeriv,
		  const double* pQuadPts,
		  const double* pQuadWts,
		  const int numDims,
		  const int numCorners,
		  const int numQuadPts);

  /** Set tolerance for minimum allowable Jacobian.
   *
   * @param tolerance Minimum allowable value for Jacobian
   */
  void jacobianTolerance(const double tolerance);

  /** Set tolerance for minimum allowable Jacobian.
   *
   * @param tolerance Minimum allowable value for Jacobian
   */
  double jacobianTolerance(void);

  /** Get number of dimensions.
   *
   * @returns Number of dimensions
   */
  int numDims(void) const;

  /** Get number of vertices in a cell.
   *
   * @returns Number of vertices in a cell
   */
  int numCorners(void) const;

  /** Get number of quadrature points.
   *
   * @param returns Number of quadrature points
   */
  int numQuadPts(void) const;

// PROTECTED MEMBERS ////////////////////////////////////////////////////
protected :

  /** Get basis functions evaluated at quadrature points.
   *
   * @returns Array of basis functions evaluated at quadrature points.
   */
  const double* basisFns(void) const;

  /** Get derivatives of basis functions evaluated at quadrature points.
   *
   * @returns Array of basis functions evaluated at quadrature points.
   */
  const double* basisFnsDeriv(void) const;

  /** Get coordinates of quadrature points.
   *
   * @returns Array of coordinates
   */
  const double* quadPts(void) const;

  /** Get weights of quadrature points.
   *
   * @returns Array of weights
   */
  const double* quadWts(void) const;

// PRIVATE MEMBERS //////////////////////////////////////////////////////
private :

  double _jacobianTol; ///< Tolernace for small Jacobian determinant
  
  /** Array of basis functions evaluated at the quadrature points.
   *
   * N1Qp1, N1Qp2, ...
   * N2Qp1, N2Qp2, ...
   *
   * size = numCorners * numQuadPts
   * index = iBasis*numQuadPts + iQuadPt
   */
  double* _pBasisFns; ///< Array of basis fns evaluated at quad pts

  /** Array of basis functions evaluated at the quadrature points.
   *
   * N1xQp1, N1yQp1, N1zQp1, N1xQp2, N1yQp2, N1zQp2, ... 
   * N2xQp1, N2yQp1, N2zQp1, N2xQp2, N2yQp2, N2zQp2, ...
   *
   * size = numCorners * numQuadPts * numDeriv
   * index = iBasis*numQuadPts*numDeriv + iQuadPt*numDeriv + iDeriv
   */
  double* _pBasisFnsDeriv;

  /** Array of coordinates of quadrature points.
   *
   * size = numQuadPts * numDims
   * index = iQuadPts*numDims + iDim
   */
  double* _pQuadPts; ///< Array of coordinates of quadrature points

  /** Array of weights of quadrature points.
   *
   * WtQp1, WtQp2, ...
   */
  double* _pQuadWts; ///< Array of weights of quadrature points

  int _numDims; ///< Number of dimensions
  int _numCorners; ///< Number of vertices in cell
  int _numQuadPts; ///< Number of quadrature points

}; // Quadrature

#include "Quadrature.icc" // inline methods

#endif // pylith_feassemble_quadrature_hh

// End of file 
