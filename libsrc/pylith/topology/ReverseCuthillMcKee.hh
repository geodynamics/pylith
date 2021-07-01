// -*- C++ -*-
//
// ======================================================================
//
// Brad T. Aagaard, U.S. Geological Survey
// Charles A. Williams, GNS Science
// Matthew G. Knepley, University at Buffalo
//
// This code was developed as part of the Computational Infrastructure
// for Geodynamics (http://geodynamics.org).
//
// Copyright (c) 2010-2021 University of California, Davis
//
// See LICENSE.md for license information.
//
// ======================================================================
//

/**
 * @file libsrc/topology/ReverseCuthillMcKee.hh
 *
 * @brief Interface to PETSc reverse Cuthill-McKee reordering.
 */

#if !defined(pylith_topology_reversecuthillmckee_hh)
#define pylith_topology_reversecuthillmckee_hh

// Include directives ---------------------------------------------------
#include "topologyfwd.hh" // forward declarations

// ReverseCuthillMcKee --------------------------------------------------
/// Interface to PETSc reverse Cuthill-McKee reordering.
class pylith::topology::ReverseCuthillMcKee
{ // ReverseCuthillMcKee

// PUBLIC MEMBERS ///////////////////////////////////////////////////////
public :

  /** Reorder vertices and cells of mesh using PETSc routines
   * implementing reverse Cuthill-McKee algorithm.
   *
   * @param mesh PyLith finite-element mesh.
   */
  static
  void reorder(topology::Mesh* mesh);

}; // ReverseCuthillMcKee

#endif // pylith_topology_reversecuthillmckee_hh


// End of file 


