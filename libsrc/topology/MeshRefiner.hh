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

#include "pylith/utils/sievetypes.hh" // USES PETSc Mesh

namespace pylith {
  namespace topology {
    class MeshRefiner;
    class TestMeshRefiner;
  } // topology

  namespace meshio {
    class DataWriter;
  } // meshio
} // pylith

namespace spatialdata {
  namespace geocoords {
    class CoordSys;
  } // geocoords
} // spatialdata

class pylith::topology::MeshRefiner
{ // MeshRefiner
  friend class TestMeshRefiner; // unit testing

// PUBLIC MEMBERS ///////////////////////////////////////////////////////
public :

  /// Constructor
  MeshRefiner(void);

  /// Destructor
  ~MeshRefiner(void);

  /** Write refined mesh.
   *
   * @param writer Data writer for partition information.
   * @param mesh Refined mesh.
   * @param cs Coordinate system for mesh.
   */
  static
  void write(meshio::DataWriter* const writer,
	     const ALE::Obj<Mesh>& mesh,
	     const spatialdata::geocoords::CoordSys* cs);

// NOT IMPLEMENTED //////////////////////////////////////////////////////
private :

  MeshRefiner(const MeshRefiner&); ///< Not implemented
  const MeshRefiner& operator=(const MeshRefiner&); ///< Not implemented

}; // MeshRefiner

#endif // pylith_topology_meshrefiner_hh


// End of file 
