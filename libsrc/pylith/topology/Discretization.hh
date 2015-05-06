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
 * @file libsrc/feassemble/Discretization.hh
 *
 * @brief C++ class for finite-element discretization of a field.
 */

#if !defined(pylith_feassemble_discretization_hh)
#define pylith_feassemble_discretization_hh

// Include directives ---------------------------------------------------
#include "pylith/topology/topologyfwd.hh" // forward declarations

#include "pylith/utils/petscfwd.h" // USES PetscFE

// Quadrature -----------------------------------------------------------
/** @brief C++ class for finite-element discretization.
 *
 * Discretization parameters collected from:
 *   + User
 *     - order of basis
 *     - continuity of basis
 *     - quadrature order
 *   + Mesh
 *     - spatial dimension
 *     - type of cell: simplex (tri/tet) or quad/hex
 *   + Field
 *     - number of components
 */
class pylith::topology::Discretization
{ // Quadrature
  friend class TestDiscretization; // unit testing

// PUBLIC METHODS ///////////////////////////////////////////////////////
public :

  /// Constructor
  Discretization(void);

  /// Destructor
  ~Discretization(void);

  /** Set order for basis functions.
   *
   * @param value Polynomial order for basis functions.
   */
  void basisOrder(const int value);

  /** Get order of basis functions.
   *
   * @returns Polynomial order of basis functions.
   */
  int basisOrder(void) const;

  /** Set basis continuity flag.
   *
   * @param value True if basis should be continuous.
   */
  void isBasisContinuous(const bool value);

  /** Get basis continuity flag.
   *
   * @returns True if basis should be continuous.
   */
  int isBasisContinuous(void) const;

  /** Set order of quadrature scheme.
   *
   * @param value Order of quadrature scheme.
   */
  void quadratureOrder(const int value);

  /** Get order of quadrature scheme.
   *
   * @returns Order of quadrature scheme.
   */
  int quadratureOrder(void) const;

  /** Create PetscFE object for discretization.
   *
   * @param dm PetscDM for finite-element mesh.
   * @param numComponents Number of components in field.
   *
   * @returns PetscFE object.
   */
  PetscFE createFE(const PetscDM dm,
		   const int numComponents);

// PRIVATE MEMBERS //////////////////////////////////////////////////////
private :

  int _basisOrder; ///< Order of 
  int _quadOrder; ///< Order of quadrature scheme.
  bool _basisContinuity; ///< Flag if basis is continuous.

// NOT IMPLEMENTED //////////////////////////////////////////////////////
private :

  Discretization(const Discretization&); ///< Not implemented.
  const Discretization& operator=(const Discretization&); ///< Not implemented

}; // Discretization

#endif // pylith_feassemble_discretization_hh


// End of file 
