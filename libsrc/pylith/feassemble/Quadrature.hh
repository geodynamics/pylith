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
// Copyright (c) 2010-2012 University of California, Davis
//
// See COPYING for license information.
//
// ======================================================================
//

/**
 * @file libsrc/feassemble/Quadrature.hh
 *
 * @brief Abstract base class for integrating over finite-elements
 * using quadrature.
 */

#if !defined(pylith_feassemble_quadrature_hh)
#define pylith_feassemble_quadrature_hh

// Include directives ---------------------------------------------------
#include "QuadratureRefCell.hh" // ISA QuadratureRefCell

#include "pylith/topology/topologyfwd.hh" // forward declarations

#include "pylith/utils/array.hh" // HASA scalar_array

// Quadrature -----------------------------------------------------------
/** @brief Abstract base class for integrating over finite-elements
 * using quadrature.
 *
 * This object contains the informatio needed to perform numerical
 * quadrature over a finite-element cell. It inherits quadrature
 * information over the reference cell from the QuadratureRefCell
 * object.
 *
 * Given a cell this object will compute the cell's Jacobian, the
 * determinant of the Jacobian, the inverse of the Jacobian, and the
 * coordinates in the domain of the cell's quadrature points. The
 * Jacobian and its inverse are computed at the quadrature points.
 */
template<typename mesh_type>
class pylith::feassemble::Quadrature : public QuadratureRefCell
{ // Quadrature
  friend class TestQuadrature; // unit testing

// PUBLIC METHODS ///////////////////////////////////////////////////////
public :

  /// Constructor
  Quadrature(void);

  /// Destructor
  ~Quadrature(void);

  /** Copy constructor.
   *
   * @param q Quadrature to copy
   */
  Quadrature(const Quadrature& q);

  /// Deallocate PETSc and local data structures.
  void deallocate(void);
  
  /// Setup quadrature engine.
  void initializeGeometry(void);

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
  const scalar_array& quadPts(void) const;

  /** Get derivatives of basis fns evaluated at quadrature points.
   *
   * @returns Array of derivatives of basis fns evaluated at
   * quadrature points
   */
  const scalar_array& basisDeriv(void) const;

  /** Get Jacobians evaluated at quadrature points.
   *
   * @returns Array of Jacobian inverses evaluated at quadrature points.
   */
  const scalar_array& jacobian(void) const;

  /** Get determinants of Jacobian evaluated at quadrature points.
   *
   * @returns Array of determinants of Jacobian evaluated at quadrature pts
   */
  const scalar_array& jacobianDet(void) const;

  /** Get precomputed geometry fields.
   *
   * @returns Geometry fields.
   */
  const topology::Fields<topology::Field<mesh_type> >*
  geometryFields(void) const;

  /** Compute geometric quantities for each cell.
   *
   * @param mesh Finite-element mesh
   * @param cells Finite-element cells for geometry.
   */
  void computeGeometry(const mesh_type& mesh,
             const ALE::Obj<typename mesh_type::SieveMesh::label_sequence>& cells);

  /** Retrieve precomputed geometric quantities for a cell.
   *
   * @param mesh Finite-element mesh
   * @param cell Finite-element cell
   */
  void retrieveGeometry(const typename mesh_type::SieveMesh::point_type& cell);

  /// Deallocate temporary storage.
  void clear(void);

  /** Compute geometric quantities for cell.
   *
   * @param coordinatesCell Coordinates of vertices in cell.
   * @param cell Finite-element cell
   */
  void computeGeometry(const scalar_array& coordinatesCell,
		       const int cell);

// PRIVATE MEMBERS //////////////////////////////////////////////////////
private :

  QuadratureEngine* _engine; ///< Quadrature geometry engine.

  /// Fields for precomputing geometry information.
  topology::Fields<topology::Field<mesh_type> >* _geometryFields;

  bool _checkConditioning; ///< True if checking for ill-conditioning.

// NOT IMPLEMENTED //////////////////////////////////////////////////////
private :

  const Quadrature& operator=(const Quadrature&); ///< Not implemented

}; // Quadrature

#include "Quadrature.icc" // inline methods
#include "Quadrature.cc" // template methods

#endif // pylith_feassemble_quadrature_hh


// End of file 
