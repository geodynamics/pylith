// -*- C++ -*-
//
// ======================================================================
//
//                           Brad T. Aagaard
//                        U.S. Geological Survey
//
// {LicenseText}
//
// ======================================================================
//

/**
 * @file pylith/topology/MeshOps.hh
 *
 * @brief Simple operations on a Mesh object.
 */

#if !defined(pylith_topology_meshops_hh)
#define pylith_topology_meshops_hh

// Include directives ---------------------------------------------------
#include "topologyfwd.hh" // forward declarations

// MeshOps --------------------------------------------------------------
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


// NOT IMPLEMENTED //////////////////////////////////////////////////////
private :

  MeshOps(void); ///< Not Implemented
  MeshOps(const MeshOps&); ///< Not implemented
  const MeshOps& operator=(const MeshOps&); ///< Not implemented


}; // MeshOps

#endif // pylith_topology_meshops_hh


// End of file 
