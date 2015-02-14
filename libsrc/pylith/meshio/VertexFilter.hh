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
 * @file libsrc/meshio/VertexFilter.hh
 *
 * @brief C++ object for filtering vertex fields when outputing
 * finite-element data.
 */

#if !defined(pylith_meshio_vertexfilter_hh)
#define pylith_meshio_vertexfilter_hh

// Include directives ---------------------------------------------------
#include "meshiofwd.hh" // forward declarations

#include "pylith/topology/topologyfwd.hh" // USES Field

// VertexFilter ---------------------------------------------------------
/** @brief C++ object for filtering vertex fields when outputing
 * finite-element data.
 */
class pylith::meshio::VertexFilter
{ // VertexFilter

// PUBLIC METHODS ///////////////////////////////////////////////////////
public :

  /// Constructor
  VertexFilter(void);

  /// Destructor
  virtual
  ~VertexFilter(void);

  /// Deallocate PETSc and local data structures.
  virtual
  void deallocate(void);
  
  /** Create copy of filter.
   *
   * @returns Copy of filter.
   */
  virtual
  VertexFilter* clone(void) const = 0;

  /** Filter field over vertices of a mesh.
   *
   * @param fieldIn Field to filter.
   */
  virtual
  topology::Field&
  filter(const topology::Field& fieldIn) = 0;

// PROTECTED METHODS ////////////////////////////////////////////////////
protected :

  /** Copy constructor.
   *
   * @param f Filter to copy.
   * @returns Pointer to this.
   */
  VertexFilter(const VertexFilter& f);

private :
  /** operator=.
  *
  * @param f Filter to copy.
  * @returns Copy of filter.
  */
  const VertexFilter& operator=(const VertexFilter& f);

}; // VertexFilter

#endif // pylith_meshio_vertexfilter_hh


// End of file 
