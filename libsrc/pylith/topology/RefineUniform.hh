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
 * @file libsrc/topology/RefineUniform.hh
 *
 * @brief Object for managing uniform global mesh refinement.
 */

#if !defined(pylith_topology_refineuniform_hh)
#define pylith_topology_refineuniform_hh

// Include directives ---------------------------------------------------
#include "topologyfwd.hh" // forward declarations

// RefineUniform --------------------------------------------------------
/// Object for managing uniform global mesh refinement.
class pylith::topology::RefineUniform
{ // RefineUniform
  friend class TestRefineUniform; // unit testing

// PUBLIC MEMBERS ///////////////////////////////////////////////////////
public :

  /// Constructor
  RefineUniform(void);

  /// Destructor
  ~RefineUniform(void);

  /** Refine mesh.
   *
   * @param newMesh Refined mesh (result).
   * @param mesh Mesh to refine.
   * @param levels Number of levels to refine.
   */
  void refine(Mesh* const newMesh,
	      const Mesh& mesh,
	      const int levels =2);

// NOT IMPLEMENTED //////////////////////////////////////////////////////
private :

  RefineUniform(const RefineUniform&); ///< Not implemented
  const RefineUniform& operator=(const RefineUniform&); ///< Not implemented

}; // RefineUniform

#endif // pylith_topology_refineuniform_hh

 
// End of file 
