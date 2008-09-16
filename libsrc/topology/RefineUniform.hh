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
 * @file pylith/topology/RefineUniform.hh
 *
 * @brief Object for managing uniform global mesh refinement.
 */

#if !defined(pylith_topology_refineuniform_hh)
#define pylith_topology_refineuniform_hh

#include "MeshRefiner.hh" // ISA MeshRefiner

namespace pylith {
  namespace topology {
    class RefineUniform;
    class TestRefineUniform;
  } // topology
} // pylith

class pylith::topology::RefineUniform : public MeshRefiner
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
  void refine(ALE::Obj<Mesh>* const newMesh,
	      const ALE::Obj<Mesh>& mesh,
	      const int levels);

// NOT IMPLEMENTED //////////////////////////////////////////////////////
private :

  RefineUniform(const RefineUniform&); ///< Not implemented
  const RefineUniform& operator=(const RefineUniform&); ///< Not implemented

}; // RefineUniform

#endif // pylith_topology_refineuniform_hh


// End of file 
