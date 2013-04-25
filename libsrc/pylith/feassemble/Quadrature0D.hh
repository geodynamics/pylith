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
// Copyright (c) 2010-2013 University of California, Davis
//
// See COPYING for license information.
//
// ======================================================================
//

/**
 * @file libsrc/feassemble/Quadrature0D.hh
 *
 * @brief Quadrature for 0-D finite-elements.
 */

#if !defined(pylith_feassemble_quadrature0d_hh)
#define pylith_feassemble_quadrature0d_hh

// Include directives ---------------------------------------------------
#include "QuadratureEngine.hh"

// Quadrature0D ---------------------------------------------------------
/** @brief Quadrature for 0-D finite-elements.
 *
 * Need Quadrature in 0-D for integration of boundary condition for
 * 1-D meshes.
 */
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
   * @param coordinatesCell Coordinates of cell's vertices.
   * @param cell Finite-element cell
   */
  void computeGeometry(const scalar_array& coordinatesCell,
		       const int cell);

  /** Compute geometric quantities for a cell at quadrature points.
   *
   * @param coordinatesCell Array of coordinates of cell's vertices.
   * @param coordinatesSize Size of coordinates array.
   * @param cell Finite-element cell
   */
  void computeGeometry(const PylithScalar* coordinatesCell,
		       const int coordinatesSize,
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
