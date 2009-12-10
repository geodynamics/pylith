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
 * @file libsrc/feassemble/QuadratureEngine.hh
 *
 * @brief Abstract base class for quadrature computation engine.
 */

#if !defined(pylith_feassemble_quadratureengine_hh)
#define pylith_feassemble_quadratureengine_hh

// Include directives ---------------------------------------------------
#include "feassemblefwd.hh" // forward declaration
#include "pylith/utils/array.hh" // USES double_array

// Quadrature0D ---------------------------------------------------------
/// Abstract base class for quadrature computation engine.
class pylith::feassemble::QuadratureEngine
{ // QuadratureEngine
  friend class TestQuadratureEngine;

// PUBLIC MEMBERS ///////////////////////////////////////////////////////
public :

  /** Constructor.
   *
   * @param q Quadrature information for reference cell.
   */
  QuadratureEngine(const QuadratureRefCell& q);

  /// Destructor
  virtual
  ~QuadratureEngine(void);

  /// Deallocate PETSc and local data structures.
  virtual
  void deallocate(void);
  
  /** Create a copy of this object.
   *
   * @returns Copy of this.
   */
  virtual
  QuadratureEngine* clone(void) const = 0;

  /** Get coordinates of quadrature points in cell (NOT reference cell).
   *
   * @returns Array of coordinates of quadrature points in cell
   */
  const double_array& quadPts(void) const;

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

  /// Allocate cell buffers.
  void initialize(void);

  /// Fill cell buffers with zeros.
  void zero(void);

  /** Compute geometric quantities for a cell at quadrature points.
   *
   * @param coordinatesCell Coordinates of cell's vertices.
   * @param cell Finite-element cell
   */
  virtual
  void computeGeometry(const double_array& coordinatesCell,
		       const int cell) = 0;

// PROTECTED METHODS ////////////////////////////////////////////////////
protected :

  /** Copy constructor.
   *
   * @param q QuadratureEngine to copy
   */
  QuadratureEngine(const QuadratureEngine& q);

  /* Check determinant of Jacobian against minimum allowable value.
   *
   * @param det Value of determinant of Jacobian
   * @param cell Label of finite-element cell
   */
  void _checkJacobianDet(const double det,
			 const int cell) const;

// PROTECTED MEMBERS ////////////////////////////////////////////////////
protected :

  /** Buffers for cell data */
  double_array _quadPts; ///< Coordinates of quad pts.
  double_array _jacobian; ///< Jacobian at quad pts;
  double_array _jacobianDet; ///< |J| at quad pts.
  double_array _jacobianInv; /// Inverse of Jacobian at quad pts.
  double_array _basisDeriv; ///< Deriv. of basis fns at quad pts.

  const QuadratureRefCell& _quadRefCell;

// NOT IMPLEMENTED //////////////////////////////////////////////////////
private :

  /// Not implemented.
  const QuadratureEngine& operator=(const QuadratureEngine&);

}; // QuadratureEngine

#include "QuadratureEngine.icc" // inline methods

#endif // pylith_feassemble_quadratureengine_hh

// End of file 
