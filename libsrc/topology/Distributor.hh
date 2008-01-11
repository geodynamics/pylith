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

#include "pylith/utils/sievefwd.hh" // USES PETSc Mesh

#include <ALE.hh>

namespace pylith {
  namespace topology {
    class Distributor;
    class TestDistributor;
  } // topology

  namespace meshio {
    class SolutionIO;
  } // meshio
} // pylith

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
  void distribute(ALE::Obj<ALE::Mesh>* const newMesh,
		  const ALE::Obj<ALE::Mesh>& mesh,
		  const char* partitioner);

  /** Write partitioning info for distributed mesh.
   *
   * @param mesh Distributed mesh.
   * @param writer Writer for partition information.
   */
  static
  void write(const ALE::Obj<ALE::Mesh>& mesh,
	     meshio::SolutionIO* const writer);

// NOT IMPLEMENTED //////////////////////////////////////////////////////
private :

  /// Not implemented
  Distributor(const Distributor&);

  /// Not implemented
  const Distributor& operator=(const Distributor&);

}; // Distributor

#endif // pylith_topology_distributor_hh


// End of file 
