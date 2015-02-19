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
 * @file libsrc/meshio/VertexFilterVecNorm.hh
 *
 * @brief C++ object for computing vector norms for fields over
 * vertices when outputing finite-element data.
 */

#if !defined(pylith_meshio_cellfiltervecnorm_hh)
#define pylith_meshio_cellfiltervecnorm_hh

// Include directives ---------------------------------------------------
#include "VertexFilter.hh" // ISA VertexFilter

// VertexFilterVecNorm --------------------------------------------------
/** @brief C++ object for computing vector norms for fields over
 * vertices when outputing finite-element data.
 */
class pylith::meshio::VertexFilterVecNorm : public VertexFilter
{ // VertexFilterVecNorm

// PUBLIC METHODS ///////////////////////////////////////////////////////
public :

  /// Constructor
  VertexFilterVecNorm(void);

  /// Destructor
  ~VertexFilterVecNorm(void);

  /** Create copy of filter.
   *
   * @returns Copy of filter.
   */
  VertexFilter* clone(void) const;

  /// Deallocate PETSc and local data structures.
  void deallocate(void);
  
  /** Filter vertex field.
   *
   * @param fieldIn Field to filter.
   */
  topology::Field&
  filter(const topology::Field& fieldIn);

// PROTECTED METHODS ////////////////////////////////////////////////////
protected :

  /** Copy constructor.
   *
   * @param f Filter to copy.
   * @returns Pointer to this.
   */
  VertexFilterVecNorm(const VertexFilterVecNorm& f);

// NOT IMPLEMENTED //////////////////////////////////////////////////////
private :

  /// Not implemented.
  const VertexFilterVecNorm& operator=(const VertexFilterVecNorm&);

// PRIVATE MEMBERS //////////////////////////////////////////////////////
private :

  topology::Field* _fieldVecNorm; ///< Filtered vertex field

}; // VertexFilterVecNorm

#endif // pylith_meshio_cellfiltervecnorm_hh


// End of file 
