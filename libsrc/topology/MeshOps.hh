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
 * @brief Temporary object for doing operations on a PETSc
 * mesh. Object will be replaced by a PyLith Mesh object that inherits
 * or templates over the PETSc mesh object.
 */

#if !defined(pylith_topology_meshops_hh)
#define pylith_topology_meshops_hh

#include "pylith/utils/sievetypes.hh" // USES PETSc Mesh

namespace pylith {
  namespace topology {
    class MeshOps;
    class TestMeshOps;
  } // topology
} // pylith

class pylith::topology::MeshOps
{ // MeshOps
  friend class TestMeshOps; // unit testing

// PUBLIC MEMBERS ///////////////////////////////////////////////////////
public :

  /** Check to make sure material id of every cell matches the id of
   *  one of the materials.
   *
   * @param mesh PETSc mesh.
   * @param materialIds Array of ids for all materials and cohesive cell interfaces.
   * @param numMaterials Size of array.
   */
  static
  void checkMaterialIds(const ALE::Obj<Mesh>& mesh,
			int* materialIds,
			const int numMaterials);


// NOT IMPLEMENTED //////////////////////////////////////////////////////
private :

  /// Not implemented
  MeshOps(const MeshOps&);

  /// Not implemented
  const MeshOps& operator=(const MeshOps&);


}; // MeshOps

#endif // pylith_topology_meshops_hh


// End of file 
