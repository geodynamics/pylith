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
 * @file pylith/topology/MeshRefiner.hh
 *
 * @brief Object for managing mesh refinement.
 */

#if !defined(pylith_topology_meshrefiner_hh)
#define pylith_topology_meshrefiner_hh

// Include directives ---------------------------------------------------
#include "topologyfwd.hh" // forward declarations

// MeshRefiner ----------------------------------------------------------
class pylith::topology::MeshRefiner
{ // MeshRefiner
  friend class TestMeshRefiner; // unit testing

// PUBLIC MEMBERS ///////////////////////////////////////////////////////
public :

  /// Constructor
  MeshRefiner(void);

  /// Destructor
  virtual
  ~MeshRefiner(void);

  /** Refine mesh.
   *
   * @param newMesh Refined mesh (result).
   * @param mesh Mesh to refine.
   * @param levels Number of levels to refine.
   * @param fields Solution fields.
   */
  virtual
  void refine(Mesh* const newMesh,
	      const Mesh& mesh,
	      const int levels =1,
	      const SolutionFields* fields =0) = 0;

// NOT IMPLEMENTED //////////////////////////////////////////////////////
private :

  MeshRefiner(const MeshRefiner&); ///< Not implemented
  const MeshRefiner& operator=(const MeshRefiner&); ///< Not implemented

}; // MeshRefiner

#endif // pylith_topology_meshrefiner_hh


// End of file 
