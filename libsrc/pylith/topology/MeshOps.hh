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
 * @file libsrc/topology/MeshOps.hh
 *
 * @brief Simple operations on a Mesh object.
 */

#if !defined(pylith_topology_meshops_hh)
#define pylith_topology_meshops_hh

// Include directives ---------------------------------------------------
#include "topologyfwd.hh" // forward declarations

// MeshOps --------------------------------------------------------------
/// Simple operations on a Mesh object.
class pylith::topology::MeshOps
{ // MeshOps
  friend class TestMeshOps; // unit testing

// PUBLIC MEMBERS ///////////////////////////////////////////////////////
public :

  /** Check to make sure material id of every cell matches the id of
   *  one of the materials.
   *
   * @param mesh Finite-element mesh.
   * @param materialIds Array of ids for all materials and cohesive
   * cell interfaces.
   * @param numMaterials Size of array.
   */
  static
  void checkMaterialIds(const Mesh& mesh,
			int* const materialIds,
			const int numMaterials);

  static
  int numMaterialCells(const Mesh& mesh, int materialId);


// NOT IMPLEMENTED //////////////////////////////////////////////////////
private :

  MeshOps(void); ///< Not Implemented
  MeshOps(const MeshOps&); ///< Not implemented
  const MeshOps& operator=(const MeshOps&); ///< Not implemented


}; // MeshOps

#endif // pylith_topology_meshops_hh


// End of file 
