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
 * This object contains the informatio needed to perform numerical
 * quadrature over a finite-element cell. It inherits quadrature
 * information over the reference cell from the QuadratureBase object.

 * Given a cell this object will compute the cell's Jacobian, the
 * determinant of the Jacobian, the inverse of the Jacobian, and the
 * coordinates in the domain of the cell's quadrature points. The
 * Jacobian and its inverse are computed at the quadrature points.
 */

#if !defined(pylith_feassemble_quadrature_hh)
#define pylith_feassemble_quadrature_hh

// Include directives ---------------------------------------------------
#include "QuadratureBase.hh" // ISA QuadratureBase

#include "pylith/topology/topologyfwd.hh" // forward declarations

#include "pylith/utils/array.hh" // HASA double_array

// Quadrature -----------------------------------------------------------
template<typename mesh_type>
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

  /** Set flag for checking ill-conditioning.
   *
   * @param flag True to check for ill-conditioning, false otherwise.
   */
  void checkConditioning(const bool flag);

  /** Get flag for checking ill-conditioning.
   *
   * @returns True if checking for ill-conditioning, false otherwise.
   */
  bool checkConditioning(void) const;

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

  /// Reset the precomputation structures.
  void resetPrecomputation(void);

  /** Precompute geometric quantities for each cell.
   *
   * @param mesh Finite-element mesh
   * @param cells Finite-element cells for geometry.
   */
  void precomputeGeometry(const mesh_type& mesh,
             const ALE::Obj<mesh_type::SieveMesh::label_sequence>& cells);

  /** Retrieve precomputed geometric quantities for a cell.
   *
   * @param mesh Finite-element mesh
   * @param cell Finite-element cell
   */
  void retrieveGeometry(const ALE::Obj<mesh_type::SieveMesh>& mesh,
                        const mesh_type::SieveMesh::point_type& cell);

// PROTECTED METHODS ////////////////////////////////////////////////////
protected :

  /** Copy constructor.
   *
   * @param q Quadrature to copy
   */
  Quadrature(const Quadrature& q);

  /* Check determinant of Jacobian against minimum allowable value.
   *
   * @param det Value of determinant of Jacobian
   * @param cell Label of finite-element cell
   */
  void _checkJacobianDet(const double det,
			 const mesh_type::SieveMesh::point_type& cell) const;

  /// Set entries in geometry arrays to zero.
  void _resetGeometry(void);

  /** Compute geometric quantities for a cell at quadrature points.
   *
   * @param coordsVertices Coordinates of vertices of cell.
   * @param spaceDim Number of coordinates per vertex.
   * @param cell Finite-element cell
   */
  virtual 
  void _computeGeometry(const double* coordsVertices,
			const int spaceDim,
			const mesh_type::SieveMesh::point_type& cell) = 0;

// PROTECTED MEMBERS ////////////////////////////////////////////////////
protected :

  /** Fields and visitors for precomputing geometry information for
   * cells associated with this quadrature.
   */
  /** :OPTIMIZATION:
   * Need to check speed of retrieving geometry using Fields and Visitors.
   * Can switch to Sieve sections and visitors if too slow.
   */
  Field<mesh_type>* _quadPtsField; ///< Coordinates of quad pts.
  Visitor<Field<mesh_type> >* _quadPtsVisitor; ///< Visitor for quad pts.

  Field<mesh_type>* _jacobianField; ///< Jacobian at quad pts.
  Visitor<Field<mesh_type> >* _jacobianVisitor; ///< Visitor for Jacobian.

  Field<mesh_type>* _jacobianDetField; ///< |J| at quad pts.
  Visitor<Field<mesh_type> >* _jacobianDetVisitor; ///< Visitor for |J|.

  Field<mesh_type>* _basisDerivField; ///< Deriv. of basis fns at quad pts.
  Visitor<Field<mesh_type> >* _basisDerivVisitor; ///< Visitor for derivatives.

  bool _checkConditioning; ///< True if checking for ill-conditioning.

// NOT IMPLEMENTED //////////////////////////////////////////////////////
private :

  const Quadrature& operator=(const Quadrature&); ///< Not implemented

}; // Quadrature

#include "Quadrature.icc" // inline methods
#include "Quadrature.cc" // template methods

#endif // pylith_feassemble_quadrature_hh


// End of file 
