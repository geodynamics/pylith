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
 * @file libsrc/feassemble/Quadrature2Din3D.hh
 *
 * @brief Quadrature for 2-D finite-elements in 3-D space.
 */

#if !defined(pylith_feassemble_quadrature2din3d_hh)
#define pylith_feassemble_quadrature2din3d_hh

// Include directives ---------------------------------------------------
#include "QuadratureEngine.hh"

// Quadrature2Din3D -----------------------------------------------------
/// Quadrature for 2-D finite-elements in 3-D space.
class pylith::feassemble::Quadrature2Din3D : public QuadratureEngine
{ // Quadrature2Din3D
  friend class TestQuadrature2Din3D; // unit testing

// PUBLIC MEMBERS ///////////////////////////////////////////////////////
public :

  /** Constructor.
   *
   * @param q Quadrature information for reference cell.
   */
  Quadrature2Din3D(const QuadratureRefCell& q);

  /// Destructor
  ~Quadrature2Din3D(void);

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
  Quadrature2Din3D(const Quadrature2Din3D& q);

// PRIVATE METHODS //////////////////////////////////////////////////////
private :

  /// Not implemented
  const Quadrature2Din3D& operator=(const Quadrature2Din3D&);

}; // Quadrature2Din3D

#include "Quadrature2Din3D.icc" // inline methods

#endif // pylith_feassemble_quadrature2din3d_hh


// End of file 
