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
// Copyright (c) 2010-2011 University of California, Davis
//
// See COPYING for license information.
//
// ======================================================================
//

/**
 * @file libsrc/feassemble/Quadrature3D.hh
 *
 * @brief Quadrature for 3-D finite-elements.
 */

#if !defined(pylith_feassemble_quadrature3d_hh)
#define pylith_feassemble_quadrature3d_hh

// Include directives ---------------------------------------------------
#include "QuadratureEngine.hh"

// Quadrature0D ---------------------------------------------------------
/// Quadrature for 3-D finite-elements in 3-D space.
class pylith::feassemble::Quadrature3D : public QuadratureEngine
{ // Quadrature3D
  friend class TestQuadrature3D; // unit testing

// PUBLIC MEMBERS ///////////////////////////////////////////////////////
public :

  /** Constructor.
   *
   * @param q Quadrature information for reference cell.
   */
  Quadrature3D(const QuadratureRefCell& q);

  /// Destructor
  ~Quadrature3D(void);

  /// Create a copy of this object.
  QuadratureEngine* clone(void) const;

  /** Compute geometric quantities for a cell at quadrature points.
   *
   * @param coordinatesCell Coordinates of cell's vertices.
   * @param cell Finite-element cell
   */
  void computeGeometry(const scalar_array& coordinatesCell,
		       const int cell);

// PROTECTED METHODS ////////////////////////////////////////////////////
protected :

  /** Copy constructor.
   *
   * @param q Quadrature to copy
   */
  Quadrature3D(const Quadrature3D& q);

// PRIVATE METHODS //////////////////////////////////////////////////////
private :

  /// Not implemented
  const Quadrature3D& operator=(const Quadrature3D&);

}; // Quadrature3D

#include "Quadrature3D.icc" // inline methods

#endif // pylith_feassemble_quadrature3d_hh


// End of file 
