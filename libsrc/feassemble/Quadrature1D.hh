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
// Copyright (c) 2010 University of California, Davis
//
// See COPYING for license information.
//
// ======================================================================
//

/**
 * @file libsrc/feassemble/Quadrature1D.hh
 *
 * @brief Quadrature for 1-D finite-elements.
 */

#if !defined(pylith_feassemble_quadrature1d_hh)
#define pylith_feassemble_quadrature1d_hh

// Include directives ---------------------------------------------------
#include "QuadratureEngine.hh" // ISA Quadrature

// Quadrature1D ---------------------------------------------------------
/// Quadrature for 1-D finite-elements.
class pylith::feassemble::Quadrature1D : public QuadratureEngine
{ // Quadrature1D
  friend class TestQuadrature1D; // unit testing

// PUBLIC MEMBERS ///////////////////////////////////////////////////////
public :

  /** Constructor.
   *
   * @param q Quadrature information for reference cell.
   */
  Quadrature1D(const QuadratureRefCell& q);

  /// Destructor
  ~Quadrature1D(void);

  /// Create a copy of this object.
  QuadratureEngine* clone(void) const;

  /** Compute geometric quantities for a cell at quadrature points.
   *
   * @param coordinatesCell Coordinates of cell's vertices.
   * @param cell Finite-element cell
   */
  void computeGeometry(const double_array& coordinatesCell,
		       const int cell);

// PRIVATE METHODS //////////////////////////////////////////////////////
private :

  /** Copy constructor.
   *
   * @param q Quadrature to copy
   */
  Quadrature1D(const Quadrature1D& q);

// NOT IMPLEMENTED //////////////////////////////////////////////////////
private :

  const Quadrature1D& operator=(const Quadrature1D&); ///< Not implemented

}; // Quadrature1D

#include "Quadrature1D.icc" // inline methods

#endif // pylith_feassemble_quadrature1d_hh


// End of file 
