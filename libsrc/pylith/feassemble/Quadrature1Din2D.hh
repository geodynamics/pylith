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
 * @file libsrc/feassemble/Quadrature1Din2D.hh
 *
 * @brief Quadrature for 1-D finite-elements in 2-D space.
 */

#if !defined(pylith_feassemble_quadrature1din2d_hh)
#define pylith_feassemble_quadrature1din2d_hh

// Include directives ---------------------------------------------------
#include "QuadratureEngine.hh"

// Quadrature1Din2D -----------------------------------------------------
/// Quadrature for 1-D finite-elements in 2-D space.
class pylith::feassemble::Quadrature1Din2D : public QuadratureEngine
{ // Quadrature1Din2D
  friend class TestQuadrature1Din2D; // unit testing

// PUBLIC MEMBERS ///////////////////////////////////////////////////////
public :

  /** Constructor.
   *
   * @param q Quadrature information for reference cell.
   */
  Quadrature1Din2D(const QuadratureRefCell& q);

  /// Destructor
  ~Quadrature1Din2D(void);

  /// Create a copy of this object.
  QuadratureEngine* clone(void) const;

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
  Quadrature1Din2D(const Quadrature1Din2D& q);

// PRIVATE METHODS //////////////////////////////////////////////////////
private :

  /// Not implemented
  const Quadrature1Din2D& operator=(const Quadrature1Din2D&);

}; // Quadrature1Din2D

#include "Quadrature1Din2D.icc" // inline methods

#endif // pylith_feassemble_quadrature1din2d_hh


// End of file 
