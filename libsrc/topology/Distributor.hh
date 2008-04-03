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
 * @file pylith/topology/Distributor.hh
 *
 * @brief Object for managing distribution of mesh among processors.
 */

#if !defined(pylith_topology_distributor_hh)
#define pylith_topology_distributor_hh

#include "pylith/utils/sievetypes.hh" // USES PETSc Mesh

namespace pylith {
  namespace topology {
    class Distributor;
    class TestDistributor;
  } // topology

  namespace meshio {
    class OutputManager;
  } // meshio
} // pylith

namespace spatialdata {
  namespace geocoords {
    class CoordSys;
  } // geocoords
} // spatialdata

class pylith::topology::Distributor
{ // Distributor
  friend class TestDistributor; // unit testing

// PUBLIC MEMBERS ///////////////////////////////////////////////////////
public :

  /// Constructor
  Distributor(void);

  /// Destructor
  ~Distributor(void);

  /** Distribute mesh among processors.
   *
   * @param newMesh Distributed mesh (result).
   * @param mesh Mesh to distribute.
   * @param partitioner Name of partitioner to use in distributing mesh.
   */
  static
  void distribute(ALE::Obj<Mesh>* const newMesh,
		  const ALE::Obj<Mesh>& mesh,
		  const char* partitioner);

  /** Write partitioning info for distributed mesh.
   *
   * @param output Output manager for partition information.
   * @param mesh Distributed mesh.
   * @param cs Coordinate system for mesh.
   */
  static
  void write(meshio::OutputManager* const output,
	     const ALE::Obj<Mesh>& mesh,
	     const spatialdata::geocoords::CoordSys* cs);

// NOT IMPLEMENTED //////////////////////////////////////////////////////
private :

  /// Not implemented
  Distributor(const Distributor&);

  /// Not implemented
  const Distributor& operator=(const Distributor&);

}; // Distributor

#endif // pylith_topology_distributor_hh


// End of file 
